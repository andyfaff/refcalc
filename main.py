from flask import Flask, render_template
from bin import slitoptimiser, utils


# If `entrypoint` is not defined in app.yaml, App Engine will look for an app
# called `app` in `main.py`.
app = Flask(__name__)


defaultd = {
    "L12": 2942.5,
    "L2S": 193,
    "LS4": 280.5,
    "LpreS1": 2166,
    "resolution": 0.033,
    "footprint": 50,
    "lambdamin": 2.8,
    "lambdamax": 18.5,
    "a1": 0.5,
    "a2": 2,
    "a3": 6,
    "a4": 1,
    "d1_a4": 1,
    "d2_a4": 1,
    "SLD1": 0,
    "SLD2": 2.07,
}


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/slits")
def slits():
    dct = defaultd.copy()
    calculate_variables(dct)
    return render_template("angulator.html", d=dct)


@app.route("/sld")
def slds():
    return "shoop"


def calculate_variables(d):
    angles = [float(d["a1"]), float(d["a2"]), float(d["a3"]), float(d["a4"])]
    footprint = float(d["footprint"])
    lambdamin = float(d["lambdamin"])
    lambdamax = float(d["lambdamax"])
    L12 = float(d["L12"])
    L2S = float(d["L2S"])
    LS4 = float(d["LS4"])
    LpreS1 = float(d["LpreS1"])
    resolution = float(d["resolution"])

    d["minqvals"] = [utils.qcalc(a, lambdamax) for a in angles]
    d["maxqvals"] = [utils.qcalc(a, lambdamin) for a in angles]

    d1, d2 = slitoptimiser.slitoptimiser(
        footprint, resolution, L12=L12, L2S=L2S, verbose=False
    )
    d["slit1"] = [d1 * angles[0], d1 * angles[1], d1 * angles[2], float(d["d1_a4"])]
    d["slit2"] = [d2 * angles[0], d2 * angles[1], d2 * angles[2], float(d["d2_a4"])]
    d["actualfootprint"] = [
        slitoptimiser.actual_footprint(w1, w2, L12, L2S, a)
        for w1, w2, a in zip(d["slit1"], d["slit2"], angles)
    ]
    d["dtheta"] = [utils.div(w1, w2, L12) for w1, w2 in zip(d["slit1"], d["slit2"])]
    d["postsampleslit"] = [
        slitoptimiser.height_of_beam_after_dx(w1, w2, L12, LS4 + L2S)
        for w1, w2 in zip(d["slit1"], d["slit2"])
    ]
    d["preS1slit"] = [
        slitoptimiser.height_of_beam_after_dx(w1, w2, L12, -LpreS1)
        for w1, w2 in zip(d["slit1"], d["slit2"])
    ]
    d["Qc"] = utils.qcrit(float(d["SLD1"]), float(d["SLD2"]))


if __name__ == "__main__":
    # This is used when running locally only. When deploying to Google App
    # Engine, a webserver process such as Gunicorn will serve the app. You
    # can configure startup instructions by adding `entrypoint` to app.yaml.
    app.run(host="127.0.0.1", port=8080, debug=True)
