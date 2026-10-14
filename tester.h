#pragma once
#include <cmath>
#include <vector>
#include "image.h"
#include <fstream>

class Tester {
public:
	Tester()
		: m_compression_factors_to_test{ 0.5, 0.75, 0.9 } {
		addAnImage({"1.ppm", "papyrus"});
	    addAnImage({"2.ppm", "red_background_image"});
		addAnImage({"3.ppm", "girl_face"});
		addAnImage({"4.ppm", "stained_glass"});
		addAnImage({ "5.ppm", "leaf" });
		addAnImage({ "6.ppm", "trees" });
		addAnImage({ "7.ppm", "water" });
		addAnImage({ "8.ppm", "city" });
		addAnImage({ "9.ppm", "cartoon cat" });
		addAnImage({ "10.ppm", "cat" });

	}

	void addAnImage(Image image) {
		m_images.emplace_back(image);
	}

	void testLowPassFilters(double step_size=0.01) {
		std::ofstream out;

		out.open("Data2to4.txt");
		int count = 0;
		for (auto& image : m_images) {
			++count;
			out << count;
			out << '\n';
			for (double exponent = 2.0; exponent <= 4.001; exponent += step_size) {

				for (auto compression_factor : m_compression_factors_to_test) {


				Image test_image = image;
				Image temp = image;
				double actual_compression_factor;
				double R = calcRFromAreaAndExponent(exponent, compression_factor, image);

				test_image.forwardFourier();
				test_image.cutOutGeneral(R, exponent, actual_compression_factor);
				test_image.reverseFourier();
					
				double SSIM = calculateSSIM(temp, test_image);

				/*if (max_SSIM < SSIM) {
					max_SSIM = SSIM;
					max_SSIM_exponent = exponent;
					max_SSIM_compression_factor = actual_compression_factor;
				}*/
				out << SSIM << '\t';
				out.flush();
				/*std::cout
					<< "Exponent: "
					<< exponent
					<< "\nIntended Compression Factor: "
					<< compression_factor
					<< "\nActual Compression Factor: "
					<< actual_compression_factor
					<< "\nSSIM: "
					<< SSIM
					<< "\n\n";*/
				}
			/*out
				<< "\n\n\n Image: "
				<< image.getName()
				<< " At Compression " 
				<< max_SSIM_compression_factor 
				<< " SSIM Max : " 
				<< max_SSIM 
				<< " with Exponent: " 
				<< max_SSIM_exponent 
				<< "\n\n\n\n";*/
				out << '\n';


			}

			out << "\n\n\n";
		}
		out.close();
	}

	/*void test() {
		double compression_factor;
		Image previous = m_images[0];
		m_images[0].writePPM("original", false);
		Image forward = m_images[0];
		m_images[0].forwardFourier("fd");
		forward.boxFilter();
		forward.forwardFourier("Boxfilteredforward");
		m_images[0].cutOutGeneral(100, 2, compression_factor, true);
		m_images[0].writePPM("cut_out", true);
		m_images[0].reverseFourier("compressed");
		m_images[0].spacialQuantisation(20);
		m_images[0].writePPM("HFcut_off", false);

		forward.cutOutGeneral(300, 2, compression_factor, false);
		forward.reverseFourier("CompressedNormal");

		Image superposition = Image{ forward, m_images[0] };

		superposition.writePPM("Superposition", false);

		std::cout << "SSIM: " << calculateSSIM(superposition, previous);
	}

	void test2() {
		double compression_factor;
		Image previous = m_images[0];
		m_images[0].writePPM("original", false);

		m_images[0].writePPM("original", false);
	}

	void test3() {
		Image start = m_images[0];
		m_images[0].localisedFourier(16, 100);
		m_images[0].writePPM("Localised", false);

		std::cout << "SSIM: " << calculateSSIM(start, m_images[0]);
	
	}

	void test4() {
		Image original = m_images[0];
		m_images[0].smartLocalisedFourier(m_images[0], 0, 0, false);
		m_images[0].writePPM("SLF", false);
		std::cout << "SSIM: " << calculateSSIM(original , m_images[0]);
	}

	void test5() {
		Image original = m_images[0];
		Image original2 = m_images[0];
		m_images[0].smartLocalisedFourier(m_images[0], 0, 0, false);
		m_images[0].writePPM("Smart", false);
		std::cout << "SSIM: " << calculateSSIM(original, m_images[0]);

		original.localisedFourier(8, 150);
		original.writePPM("ReverseFourier", false);

		std::cout << "SSIM: " << calculateSSIM(original, original2);
	}*/

	void testSubdividing1() {
		std::ofstream file{ "output/_results.txt" };

		int image_num = 0;

		for (auto& image : m_images) {
			++image_num;
			std::cout << image_num << '\n';

			double cr_smart = 0;
			double cr_dumb = 0;

			Image dumb = image;
			Image original = image;

			dumb.localisedFourier(8, 150, cr_dumb, 0.9);

			double sobel = 60;
			double michelson = 0.3;

			int trial_count = 0;
			do {
				++trial_count;
				image.smartLocalisedFourier(image, 0, 0, cr_smart, sobel, michelson, false);

				if (std::abs(cr_smart - cr_dumb) >= 0.01) {
					image.setData(original);
					if (cr_smart > cr_dumb) {
						--sobel;
					}
					else {
						++sobel;
					}
				}
				std::cout << "Trial: " << trial_count << " CR Smart: " << cr_smart << " Cr Dumb: " << cr_dumb << '\n';
			} 
			while (std::abs(cr_smart - cr_dumb) >= 0.01);
			std::cout << "\n\n";

			image.writePPM("Smart", false);
			dumb.writePPM("ReverseFourier", false);

			file << cr_smart << '\t' << calculateSSIM(original, image) << '\t';
			file << cr_dumb << '\t' << calculateSSIM(dumb, original);
			file << '\n';
		}
	}
	void testSubdividing2() {
		std::ofstream file{ "output/_results.txt" };

		int image_num = 0;

		for (auto& image : m_images) {
			++image_num;
			std::cout << image_num << '\n';

			double cr_smart = 0;
			double cr_dumb = 0;

			Image dumb = image;
			Image original = image;

			dumb.localisedFourier(8, 150, cr_dumb, 0.9);

			double sobel = 60;
			double michelson = 0.5;

			image.smartLocalisedFourier(image, 0, 0, cr_smart, sobel, michelson, false);

			image.writePPM("Smart", false);
			dumb.writePPM("ReverseFourier", false);

			file << cr_smart << '\t' << calculateSSIM(original, image) << '\t';
			file << cr_dumb << '\t' << calculateSSIM(dumb, original);
			file << '\n';
		}
	}
	void testSubdividing3() {

		std::vector<double> test_crs{ 0.5, 0.75, 0.9, 0.95 };
		
		for (int j = 0; j < test_crs.size(); ++j) {
			int image_num = 0;
			std::ofstream file{ "output2/_" + std::to_string(test_crs[j]) + "results.txt"};
			for (auto& image : m_images) {

				file << m_images[image_num].getName() << "\n\n";
				for (int i = 0; i < 63; ++i) {
					std::cout << "hi";
					double cr = 0;


					Image modify = image;
					modify.localisedFourier(i + 2, 150, cr, test_crs[j]);

					if (i == 12) {
						modify.writePPM(std::to_string(i + 2), false);
					}

					file << i + 2 << '\t' << cr << '\t' << calculateSSIM(modify, image) << '\t';
					file << '\n';
					file.flush();
				}

				++image_num;
			}
		}
		
	}
	void testSubdividing4() {
		std::ofstream file{ "output/_SmartJPEGresults2.txt" };

		std::vector<double> test_crs{ 0.5, 0.75, 0.9, 0.95 };

		for (auto test_cr : test_crs) {
			for (auto& image : m_images) {
				double cr_smart = 0;

				Image smart = image;
				Image smart2 = image;

				double sobel = 80;
				double michelson = 0.6;

				smart.smartLocalisedFourier(image, 0, 0, cr_smart, sobel, michelson, test_cr, false);


				file << image.getName() << '\t';
				file << cr_smart << '\t' << calculateSSIM(smart, image) << '\t';
				file << '\n';
				file.flush();
			}
		}
	}

	void getOriginalImages() {
		for (auto& image : m_images) {
			image.writePPM("OriginalBAW", false);
		}
	}

private:
	double calcRFromAreaAndExponent(double exponent, double compression_factor, const Image& image) {
		double ellipse_scale_factor = image.getRows() / static_cast<double>(image.getCols());
		double final_area = (1.0 - compression_factor) * image.getPixelCount();

		double R = std::sqrt(final_area * std::tgamma(1 + 2 / exponent) / (4 * ellipse_scale_factor * std::pow(std::tgamma(1 + 1 / exponent), 2)));

		return R;
	}

	std::vector<Image> m_images;
	std::vector<double> m_compression_factors_to_test;
};

