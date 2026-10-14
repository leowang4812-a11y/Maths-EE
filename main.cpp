#include <iostream>
#include "kiss_fftnd.h"
#include "image.h"
#include "tester.h"

void transformationPackageDoubleFourier(Image img) {
	Image before_transform = img;
	img.writePPM("start", false);
	img.forwardFourier("frequency");
	//img.averagingOutSquares(8, img.getCols() / 7);
	//img.cutOutEllipse(img.getCols() / 7);
	img.doubleFourier(16, img.getCols() / 5, 2);
	img.writePPM("frequency_cutout", true);
	img.reverseFourier("compressed");

	std::cout << "SSIM = " << calculateSSIM(img, before_transform) << '\n';
}

void transformationPackageEllipse(Image img) {
	Image before_transform = img;
	img.writePPM("start", false);
	img.forwardFourier("frequency");
	//img.cutOutGeneral(img.getCols() / 20, 1.5);
	img.writePPM("frequency_cutout", true);
	img.reverseFourier("compressed");


	std::cout << "SSIM = " << calculateSSIM(img, before_transform) << '\n';
}

int main() {
	Tester tester;
	tester.testSubdividing4();

	return 0;
}
