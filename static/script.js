const length = document.getElementById('length');
const width = document.getElementById('width');
const nSLD1 = document.getElementById('SLD1');
const nSLD2 = document.getElementById('SLD2');
const a1 = document.getElementById('a1');
const a2 = document.getElementById('a2');
const a3 = document.getElementById('a3');
const a4 = document.getElementById('a4');
const lambdamin = document.getElementById('lambdamin');
const lambdamax = document.getElementById('lambdamax');


const hypoteneuse = document.getElementById('hypoteneuse');
const Qc = document.getElementById('Qc');
const minQa1 = document.getElementById('minQa1');
const minQa2 = document.getElementById('minQa2');
const minQa3 = document.getElementById('minQa3');
const minQa4 = document.getElementById('minQa4');
const maxQa1 = document.getElementById('maxQa1');
const maxQa2 = document.getElementById('maxQa2');
const maxQa3 = document.getElementById('maxQa3');
const maxQa4 = document.getElementById('maxQa4');

a1.addEventListener('input', (event) => {
    minQa1.textContent = q(a1.value, lambdamax.value);
    maxQa1.textContent = q(a1.value, lambdamin.value);
 });

a2.addEventListener('input', (event) => {
    minQa2.textContent = q(a2.value, lambdamax.value);
    maxQa2.textContent = q(a2.value, lambdamin.value);
});

a3.addEventListener('input', (event) => {
    minQa3.textContent = q(a3.value, lambdamax.value);
    maxQa3.textContent = q(a3.value, lambdamin.value);
});

a4.addEventListener('input', (event) => {
    minQa4.textContent = q(a4.value, lambdamax.value);
    maxQa4.textContent = q(a4.value, lambdamin.value);
});

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
    minQa1.textContent = q(a1.value, lambdamax.value);
    maxQa1.textContent = q(a1.value, lambdamin.value);
    minQa2.textContent = q(a2.value, lambdamax.value);
    maxQa2.textContent = q(a2.value, lambdamin.value);
    minQa3.textContent = q(a3.value, lambdamax.value);
    maxQa3.textContent = q(a3.value, lambdamin.value);
    minQa4.textContent = q(a4.value, lambdamax.value);
    maxQa4.textContent = q(a4.value, lambdamin.value);
};

function hypo(len, wid){
    return Math.sqrt(len**2 + wid ** 2).toFixed(1)
};

function qcrit(SLD1, SLD2){
    return Math.sqrt(16.0 * Math.PI * (SLD2 - SLD1) * 1.0e-6).toFixed(5);
};

function q(angle, wavelength){
    return (4 * Math.PI * Math.sin(angle * Math.PI / 180.) / wavelength).toFixed(4);
};