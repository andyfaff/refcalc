import tomllib
from flask import Flask, render_template, request
from bin import slitoptimiser, utils

import periodictable as pt


# If `entrypoint` is not defined in app.yaml, App Engine will look for an app
# called `app` in `main.py`.
app = Flask(__name__)


defaultd = {
    "instrument": "Platypus",
    # "L12": 2990.9,
    # "L2S": 144,
    # "LS4": 444.5,
    # "LpreS1": 2166,
    "resolution": 0.033,
    "footprint": 50,
    "lambdamin": 2.8,
    "lambdamax": 18.5,
    "a1": 0.8,
    "a2": 3.5,
    "a3": 6,
    "a4": 1,
    "d1_a4": 1,
    "d2_a4": 1,
    "SLD1": 0,
    "SLD2": 2.07,
}
with open("bin/config.toml", "rb") as f:
    instrument_config = tomllib.load(f)
defaultd.update(instrument_config["Platypus"])


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/slits", methods=["POST", "GET"])
def slits():
    dct = defaultd.copy()

    if request.method == "POST":
        dct = {k: v for k, v in request.form.items()}
        instrument = dct["instrument"]
        if instrument != defaultd["instrument"]:
            # we're changing the instrument type, update distances
            defaultd["instrument"] = instrument
            settings = instrument_config[instrument]
            dct.update(settings)
            #defaultd.update(settings)

    calculate_variables(dct)
    return render_template("angulator.html", d=dct)


@app.route("/singleslit", methods=["POST", "GET"])
def singleslit():

    if request.method == "GET":
        return f"send a post request setting variables in ', {defaultd.keys()}"

    elif request.method == "POST":

        dct = defaultd.copy()
        _form = {k: float(v) for k, v in request.form.items() if k in dct}

        dct.update(_form)

        d1, d2 = slitoptimiser.slitoptimiser(
            dct["footprint"],
            dct["resolution"],
            angle=dct["a1"],
            L12=dct["L12"],
            L2S=dct["L2S"],
            verbose=False,
        )

        postslit = slitoptimiser.height_of_beam_after_dx(
            d1, d2, dct["L12"], dct["LS4"] + dct["L2S"]
        )
        preslit = slitoptimiser.height_of_beam_after_dx(
            d1, d2, dct["L12"], -dct["LpreS1"]
        )

        return f"{(preslit[1], d1, d2, postslit[1])}"


@app.route("/sld", methods=["POST", "GET"])
def slds():
    if request.method == "GET":
        d = {
            "formula": "",
            "density": 1.0,
            "volume": 30.0,
            "volumetype": "density",
            "density_checked": "checked",
            "volume_checked": "",
            "xray_energy": 8.048,
            "neutron_wavelength": 1.8,
            "neutron_sld": 0.0 + 0.0 * 1j,
            "xray_sld": 0.0 + 0.0 * 1j,
        }
        return render_template("sldcalculator.html", d=d)

    elif request.method == "POST":
        dct = {k: v for k, v in request.form.items()}

        volumetype = dct["volumetype"]
        if volumetype == "density":
            dct["density_checked"] = "checked"
            dct["volume_checked"] = ""
        elif volumetype == "volume":
            dct["density_checked"] = ""
            dct["volume_checked"] = "checked"

        dct["volume"] = volume = float(dct["volume"])
        dct["density"] = density = float(dct["density"])

        try:
            formula = pt.formula(
                dct["formula"],
            )
        except:
            return render_template("sldcalculator.html", d=dct)

        if volumetype == "volume":
            density = formula.molecular_mass / formula.volume(a=volume, b=1, c=1)
            dct["density"] = density
        elif volumetype == "density":
            volume = formula.mass / density / pt.constants.avogadro_number * 1e24
            dct["volume"] = volume

        try:
            real, imag, mu = pt.neutron_sld(
                formula, density=density, wavelength=float(dct["neutron_wavelength"])
            )
            dct["neutron_sld"] = real + imag * 1j

            real, imag = pt.xray_sld(
                formula, density=density, energy=float(dct["xray_energy"])
            )

            dct["xray_sld"] = real + imag * 1j
        except (TypeError, AssertionError):
            pass

        return render_template("sldcalculator.html", d=dct)


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
    app.run(host="0.0.0.0", port=8080, debug=True)
