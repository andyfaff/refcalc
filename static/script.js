const length = document.getElementById('length');
const width = document.getElementById('width');
const hypoteneuse = document.getElementById('hypoteneuse');

length.addEventListener('input', (event) => {
    hypoteneuse.textContent = Math.sqrt(length.value**2 + width.value ** 2);
});

width.addEventListener('input', (event) => {
    hypoteneuse.textContent =  Math.sqrt(length.value**2 + width.value ** 2);
});

