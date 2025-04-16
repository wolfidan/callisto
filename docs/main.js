// Add your dates here (or generate them dynamically if hosted by a static generator)
const imageDates = [
  "20240401",
  "20240402",
  "20240403",
  // etc.
];

let index = imageDates.length - 1;  // Start with most recent

function updateImage() {
  const date = imageDates[index];
  const img = document.getElementById("plot");
  img.src = `images/callisto_meteoswiss_${date}.png`;
}

function prevImage() {
  if (index > 0) {
    index--;
    updateImage();
  }
}

function nextImage() {
  if (index < imageDates.length - 1) {
    index++;
    updateImage();
  }
}

document.addEventListener("DOMContentLoaded", updateImage);

