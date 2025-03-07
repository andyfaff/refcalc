const length = document.getElementById('length');
const width = document.getElementById('width');
const nSLD1 = document.getElementById('SLD1');
const nSLD2 = document.getElementById('SLD2');

const hypoteneuse = document.getElementById('hypoteneuse');
const Qc = document.getElementById('Qc');

length.addEventListener('input', (event) => {
    hypoteneuse.textContent =  hypo(length.value, width.value);
});

width.addEventListener('input', (event) => {
    hypoteneuse.textContent =  hypo(length.value, width.value);
});

nSLD1.addEventListener('input', (event) => {
    Qc.textContent = qcrit(nSLD1.value, nSLD2.value);
});

nSLD2.addEventListener('input', (event) => {
    Qc.textContent = qcrit(nSLD1.value, nSLD2.value);
});

window.onload = function() {
    hypoteneuse.textContent = hypo(length.value, width.value);
    Qc.textContent = qcrit(nSLD1.value, nSLD2.value);
};

function hypo(len, wid){
    return Math.sqrt(len**2 + wid ** 2).toFixed(1)
};

function qcrit(SLD1, SLD2){
    return Math.sqrt(16.0 * Math.PI * (SLD2 - SLD1) * 1.0e-6).toFixed(5);
};
