#pragma once
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string_view>
#include "kiss_fftnd.h"
#include <random>
#include <map>

class Image {
public:
	Image(std::string path, std::string name)
		: m_path{ path }
		, m_name{ name } {

		std::ifstream ppm{ m_path };
		std::stringstream buffer;

		buffer << ppm.rdbuf();

		std::string word;

		buffer >> word; buffer >> word;
		m_cols = std::stoi(word);
		buffer >> word;
		m_rows = std::stoi(word);
		buffer >> word;

		std::vector<kiss_fft_cpx> data{ static_cast<unsigned int>(m_rows * m_cols) };

		int index = 0, count = 0;
		std::array<int, 3> colours;

		while (buffer >> word) {
			int value = std::stoi(word);
			int count_mod_3 = count % 3;
			
			colours[count_mod_3] = value;
			if (count % 3 == 2) {
				data[index].r = 0.299 * colours[0] + 0.587 * colours[1] + 0.114 * colours[2];
				data[index].i = 0;
				++index;
			}
			++count;
		}

		if (count / 3 != m_rows * m_cols) {
			std::cerr << "Mismatch of Image Size! image_size: " << " and cols * rows = " << m_cols * m_rows;
		}
		m_data = data;
		ppm.close();
	}
	Image(Image& image1, Image& image2)
		: m_name{ image1.m_name + "_and_" + image2.m_name }
		, m_cols{ image1.m_cols }
		, m_rows{ image1.m_rows } {
		if (image1.m_data.size() != image2.m_data.size()) {
			std::cerr << "Error: Image Size Mismatch\n";
		}
		std::vector<kiss_fft_cpx> data{ static_cast<unsigned int>(m_rows * m_cols) };
		for (int i = 0; i < image1.m_data.size(); ++i) {
			data[i].r = (image1.m_data[i].r + image2.m_data[i].r);
			data[i].i = (image1.m_data[i].i + image2.m_data[i].i);
		}
		m_data = data;
	}
	Image(int cols, int rows)
		: m_cols{ cols }
		, m_rows{ rows }
		, m_name{ "square" } {
		std::vector<kiss_fft_cpx> data{ static_cast<unsigned int>(m_rows * m_cols) };
		m_data = data;
	}
	~Image() = default;

	void writePPM(std::string name_extension, bool logarithmic_scaling) {
		std::ofstream ppm{ "BAW/" + m_name + '_' + name_extension + ".ppm"};
		ppm << "P3\n" << m_cols << ' ' << m_rows << "\n255\n";
		if (logarithmic_scaling)
			for (int i = 0; i < m_data.size(); ++i) {
				int value = std::abs(m_data[i].r) + std::abs(m_data[i].i);
				value = 25 * std::log10(1 + value);
				ppm << value << ' ' << value << ' ' << value << '\n';
			}
		else
			for (int i = 0; i < m_data.size(); ++i) {
				int value = m_data[i].r;
				ppm << value << ' ' << value << ' ' << value << '\n';
			}
		ppm.close();
	}

	void forwardFourier(std::string name_extension="none") {
		int dims[2] = { m_rows, m_cols };

		kiss_fftnd_cfg cfg = kiss_fftnd_alloc(dims, 2, 0, nullptr, nullptr);
		std::vector<kiss_fft_cpx> out(m_rows * m_cols);
		kiss_fftnd(cfg, m_data.data(), out.data());

		m_data = out;

		FFTShift();
		delete cfg;
		if (name_extension == "none") return;
		else {
			writePPM(name_extension, true);
		}

	}

	void reverseFourier(std::string name_extension = "none") {
		FFTShift();
		int dims[2] = { m_rows, m_cols };

		kiss_fftnd_cfg cfg = kiss_fftnd_alloc(dims, 2, 1, nullptr, nullptr);
		std::vector<kiss_fft_cpx> out(m_rows * m_cols);
		kiss_fftnd(cfg, m_data.data(), out.data());

		int total_size = m_cols * m_rows;
		for (int i = 0; i < m_data.size(); ++i) {
			m_data[i].r = out[i].r / total_size;
			m_data[i].i = out[i].i / total_size;
		}

		delete cfg;
		if (name_extension == "none") return;
		else {
			writePPM(name_extension, false);
		}
	}
	
	void cutOutCircle(int radius) {
		int row = 0;
		int col = 0;
		for (int i = 0; i < m_data.size(); ++i) {
			if ((row-(m_rows/2))*(row-(m_rows/2)) + (col - (m_cols/2)) * (col-(m_cols/2)) > radius * radius) {
				m_data[i].r = 0;
				m_data[i].i = 0;
			}

			++col;
			if (col == m_cols) {
				++row;
				col = 0;
			}
		}
	}

	void cutOutEllipse(int radius, bool discard=false) {
		int row = 0;
		int col = 0;
		for (int i = 0; i < m_data.size(); ++i) {
			if ((row - (m_rows / 2.0)) * (row - (m_rows / 2.0)) / (static_cast<float>(m_rows) / m_cols) + (col - (m_cols / 2.0)) * (col - (m_cols / 2.0)) > radius * radius) {
				m_data[i].r = 0;
				m_data[i].i = 0;
			}
			if (discard) {
				if ((row - (m_rows / 2.0)) * (row - (m_rows / 2.0)) / (static_cast<float>(m_rows) / m_cols) + (col - (m_cols / 2.0)) * (col - (m_cols / 2.0)) < radius * radius) {
					m_data[i].r = 0;
					m_data[i].i = 0;
				}
			}

			++col;
			if (col == m_cols) {
				++row;
				col = 0;
			}
		}
	}
	//row and col_cut are inclusive so they don't cut the row / column they are set at.
	void cutOutHyperbola(double a) {
		int row = 0;
		int col = 0;
		for (int i = 0; i < m_data.size(); ++i) {
			if (std::abs((row - (m_rows / 2.0)) * (col - (m_cols / 2.0))) > a * a / 2) {
				m_data[i].r = 0;
				m_data[i].i = 0;
			}

			++col;
			if (col == m_cols) {
				++row;
				col = 0;
			}
		}
	}

	void cutOutGeneral(int radius, double exponent, double& get_compression_factor, bool keep_high_fd = false) {
		int row = 0;
		int col = 0;

		double ellipse_scale_factor = m_rows / static_cast<double>(m_cols);
		int area_remaining = 0;

		for (int i = 0; i < m_data.size(); ++i) {
			double x_component = std::pow(std::abs(row - (m_rows / 2.0)) / ellipse_scale_factor, exponent);
			double y_component = std::pow(std::abs(col - (m_cols / 2.0)), exponent);
			if (keep_high_fd) {
				if (x_component + y_component < std::pow(radius, exponent)) {
					m_data[i].r = 0;
					m_data[i].i = 0;
				}
				else {
					++area_remaining;
				}
			}
			else {
				if (x_component + y_component > std::pow(radius, exponent)) {
					m_data[i].r = 0;
					m_data[i].i = 0;
				}
				else {
					++area_remaining;
				}
			}
			++col;
			if (col == m_cols) {
				++row;
				col = 0;
			}
		}

		get_compression_factor = 1.0 - area_remaining / static_cast<double>(m_data.size());

	}

	void removeOuterEdges(int row_cut, int col_cut) {
		int row = 0;
		int col = 0;
		for (int i = 0; i < m_data.size(); ++i) {
			if (row > row_cut && col > col_cut && row < m_rows - row_cut && col < m_cols - col_cut) {
				m_data[i].r = 0;
				m_data[i].i = 0;
			}

			++col;
			if (col == m_cols) {
				++row;
				col = 0;
			}
		}
	}

	void spacialQuantisation(int cut_off) {
		for (auto& value : m_data) {
			if (value.r < cut_off) {
				value.r = 0;
				value.i = 0;
			}
		}
	}

	double calcRunLengthCompression() {
		int unique_data_points = 0;
		int prev_data_point = 256;
		for (auto& data : m_data) {
			if (data.r != prev_data_point) {
				++unique_data_points;
				prev_data_point = data.r;
			}
		}
		return unique_data_points * 2;
	}

	void boxFilter() {
		Image filtered{ m_cols, m_rows };
		for (int i = 0; i < m_data.size(); ++i) {
			int x = i % m_rows;
			int y = i / m_rows;

			int count = 0;
			double sum = 0;
			for (int a = 0; a < 5; ++a) {
				for (int b = 0; b < 5; ++b) {
					int xk = x + a;
					int yk = y + b;

					if (0 <= xk < m_cols && 0 <= yk < m_rows) {
						sum += m_data[yk * m_cols + xk].r;
						++count;
					}
				}
			}
			filtered.m_data[y * m_cols + x].r = sum / count;
		}
		m_data = filtered.m_data;
	}

	void averagingOutSquares(int side_length, double radius) {
		int square_rows = std::ceil((m_rows + 1) / static_cast<double>(side_length));
		int square_cols = std::ceil((m_cols + 1) / static_cast<double>(side_length));

		for (int i = 0; i < square_rows; ++i) {
			for (int j = 0; j < square_cols; ++j) {
				if (std::pow(i - (square_rows / 2.0), 2) / (static_cast<float>(m_rows) / m_cols) + std::pow(j - (square_cols / 2.0), 2) < std::pow(radius / side_length, 2)) {
					continue;
				}
				int corner_x = j * side_length;
				int corner_y = i * side_length;

				double mean_r = 0;
				double mean_i = 0;

				int length_x = side_length;
				int length_y = side_length;

				if (j == square_cols - 1) {
					length_x = m_cols % side_length;
				}

				if (i == square_rows - 1) {
					length_y = m_rows % side_length;
				}
					
				for (int b = 0; b < length_y; ++b) {
					for (int a = 0; a < length_x; ++a) {
						int index = (corner_y + b) * m_cols + (corner_x + a);
						mean_r += m_data[index].r;
						mean_i += m_data[index].i;

					}
				}

				mean_r /= std::pow(side_length, 2);
				mean_i /= std::pow(side_length, 2);

				for (int b = 0; b < length_y; ++b) {
					for (int a = 0; a < length_x; ++a) {
						int index = (corner_y + b) * m_cols + (corner_x + a);
						m_data[index].r = mean_r;
						m_data[index].i = mean_i;
					}
				}
			}
		}

	}

	void quantisationMatrix(int cutoff_ratio, int& remaining) {
		int indexX =  m_cols / 2;
		int indexY = m_rows / 2;
		std::vector<kiss_fft_cpx> centres{ {
			m_data[indexY * m_cols + indexX],
			m_data[indexY * m_cols + indexX - 1],
			m_data[(indexY - 1) * m_cols + indexX],
			m_data[(indexY - 1) * m_cols + indexX - 1]
		} };

		double max_value = 0;

		for (int i = 0; i < 4; ++i) {
			int value = std::sqrt(std::pow(centres[i].r, 2) + std::pow(centres[i].i, 2));
			if (value > max_value) {
				max_value = value;
			}
		}
		for (int i = 0; i < m_data.size(); ++i) {
			int value = std::sqrt(std::pow(m_data[i].r, 2) + std::pow(m_data[i].i, 2));
			if (value > max_value) {
				max_value = value;
			}
		}

		int count = 0;
		for (int i = 0; i < m_data.size(); ++i) {
			if (cutoff_ratio * (std::sqrt(std::pow(m_data[i].r, 2) + std::pow(m_data[i].i, 2))) / max_value < 1.0) {
				m_data[i].r = 0;
				m_data[i].i = 0;

			}
			else {
				++remaining;
			}
		}
	}

	void doubleFourier(int side_length, double radius, int cutoff_ratio) {
		int square_rows = std::ceil((m_rows + 1) / static_cast<double>(side_length));
		int square_cols = std::ceil((m_cols + 1) / static_cast<double>(side_length));

		for (int i = 0; i < square_rows; ++i) {
			for (int j = 0; j < square_cols; ++j) {
				if (std::pow(i - (square_rows / 2.0 - 1) , 2) / (static_cast<float>(m_rows) / m_cols) + std::pow(j - (square_cols / 2.0 - 1), 2) < std::pow(radius / side_length, 2)) {
					continue;
				}
				int corner_x = j * side_length;
				int corner_y = i * side_length;

				int length_x = side_length;
				int length_y = side_length;

				if (j == square_cols - 1) {
					length_x = m_cols % side_length;
				}

				if (i == square_rows - 1) {
					length_y = m_rows % side_length;
				}

				Image square{ length_x, length_y };
				int max_value = 0;
				for (int b = 0; b < length_y; ++b) {
					for (int a = 0; a < length_x; ++a) {
						int index = (corner_y + b) * m_cols + (corner_x + a);
						square.m_data[b * length_x + a].r = m_data[index].r;
						square.m_data[b * length_x + a].i = m_data[index].i;
					}
				}
				if (length_x > 0 && length_y > 0) {
					square.forwardFourier();
					int remaining;
					square.quantisationMatrix(cutoff_ratio, remaining);
					square.reverseFourier();
				}

				for (int b = 0; b < length_y; ++b) {
					for (int a = 0; a < length_x; ++a) {
						int index = (corner_y + b) * m_cols + (corner_x + a);
						m_data[index] = square.m_data[b * square.m_cols + a];
					}
				}
			}
		}
	}

	void localisedFourier(int side_length, int cutoff_ratio, double& compression_ratio, double target_cr) {
		int square_rows = std::ceil((m_rows + 1) / static_cast<double>(side_length));
		int square_cols = std::ceil((m_cols + 1) / static_cast<double>(side_length));

		int remaining = 0;

		for (int i = 0; i < square_rows; ++i) {
			for (int j = 0; j < square_cols; ++j) {
				int corner_x = j * side_length;
				int corner_y = i * side_length;

				int length_x = side_length;
				int length_y = side_length;

				if (j == square_cols - 1) {
					length_x = m_cols % side_length;
				}

				if (i == square_rows - 1) {
					length_y = m_rows % side_length;
				}

				Image square{ length_x, length_y };
				int max_value = 0;
				for (int b = 0; b < length_y; ++b) {
					for (int a = 0; a < length_x; ++a) {
						int index = (corner_y + b) * m_cols + (corner_x + a);
						square.m_data[b * length_x + a].r = m_data[index].r;
						square.m_data[b * length_x + a].i = m_data[index].i;
					}
				}
				if (length_x > 0 && length_y > 0) {
					square.forwardFourier();

					square.varyingQM(remaining, target_cr);
					square.reverseFourier();
				}

				for (int b = 0; b < length_y; ++b) {
					for (int a = 0; a < length_x; ++a) {
						int index = (corner_y + b) * m_cols + (corner_x + a);
						m_data[index] = square.m_data[b * square.m_cols + a];
					}
				}
			}
		}
		compression_ratio = static_cast<double>(remaining) / m_data.size();
	}

	static double calcFourierEnergy(Image& img) {
		double cf;
		img.forwardFourier("test");
		Image forward = img;
		img.cutOutGeneral(img.getRows() / 5, 2, cf, true);
		double total_high = 0;
		double total = 0;

		for (int i = 0; i < img.m_data.size(); ++i) {
			total_high += img.m_data[i].r;
			total += forward.m_data[i].r;
		}

		return total_high / total;
	}

	static double michelsonContrast(const Image& img) {
		double max = 0;
		double min = INFINITY;
		
		for (const auto& data : img.m_data) {
			if (data.r > max) {
				max = data.r;
			}
			if (data.r < min) {
				min = data.r;
			}
		}
		return (max - min) / (max + min);
	}

	static double gradientSobelCalc(const Image& img) {
		Image gradientx{ img.m_cols, img.m_rows };
		Image gradienty{ img.m_cols, img.m_rows };
		Image gradientFinal{ img.m_cols, img.m_rows };

		int mean = img.getMean();

		for (int i = 0; i < img.m_data.size(); ++i) {
			int x = i % img.m_rows;
			int y = i / img.m_rows;

			auto addToGradient = [&](int multiply, int x_trans, int y_trans, Image& gradient) {
				int index = (y + y_trans) * img.m_cols + x + x_trans;
				if (y + y_trans < 0 || y + y_trans >= img.m_rows || x + x_trans < 0 || x + x_trans >= img.m_cols) {
					gradient.m_data[i].r += mean * multiply;
				}
				else {
					gradient.m_data[i].r += img.m_data[index].r * multiply;
				}
			};

			addToGradient(-1, -1, -1, gradientx);
			addToGradient(-2, -1, 0, gradientx);
			addToGradient(-1, -1, 1, gradientx);
			addToGradient(1, 1, -1, gradientx);
			addToGradient(2, 1, 0, gradientx);
			addToGradient(1, 1, 1, gradientx);

			addToGradient(-1, -1, -1, gradienty);
			addToGradient(-2, 0, -1, gradienty);
			addToGradient(-1, 1, -1, gradienty);
			addToGradient(1, -1, 1, gradienty);
			addToGradient(2, 0, 1, gradienty);
			addToGradient(1, 1, 1, gradienty);

			gradientFinal.m_data[i].r = std::pow(gradientx.m_data[i].r, 2) + std::pow(gradienty.m_data[i].r, 2);
		}

		return std::sqrt(gradientFinal.getSD(gradientFinal.getMean()));

	}

	double runLengthEncoding(const std::vector<std::vector<int>>& after_Qmatrix) {
		std::vector infoZigZag{ m_data.size() };
		int x = 0;
		int y = 0;
		int xmax = 0;
		for (int i = 0; i < m_cols; ++i) {
			infoZigZag.emplace_back();
			if (x == 0) {
				++xmax;
				x = xmax;
				y = 0;
			}
			--x;
			++y;
		}
	}

	void varyingQM(int& count, double target_cr) {
		std::vector<std::vector<double>> Qmatrix(m_rows, std::vector<double>(m_cols));
		std::vector<std::pair<double, double>> values_after_division{};
		static double remove_carry_on = 0.0;
		int temp_count = 0;
		for (int i = 0; i < m_rows; ++i) {
			for (int j = 0; j < m_cols; ++j) {
				double Qvalue = QM16[static_cast<int>(i * 16.0 / m_rows)][static_cast<int>(j * 16.0 / m_cols)] * 4 / (256.0 / m_data.size());

				double value_after_division = (std::sqrt(std::pow(m_data[i * m_cols + j].r, 2) + std::pow(m_data[i * m_cols + j].i, 2))) / Qvalue;
				values_after_division.emplace_back(std::pair{ i * m_cols + j, value_after_division });



				//std::cout << std::sqrt(std::pow(m_data[i * m_cols + j].r, 2) + std::pow(m_data[i * m_cols + j].i, 2)) << ' ';
			}
			//std::cout << "\n";
		}
		//std::cout << "\n\n\n";
		std::sort(values_after_division.begin(), values_after_division.end(), [](auto& left, auto& right) {
			return left.second < right.second;
			});
	
		double upper_bound = m_data.size() * target_cr;

		while (remove_carry_on >= 1.0) {
			remove_carry_on -= 1;
			--upper_bound;
		}

		for (int n = 0; n < upper_bound; ++n) {

			m_data[values_after_division[n].first].r = 0;
			m_data[values_after_division[n].first].i = 0;

			++temp_count;
		}

		remove_carry_on += ((double)(temp_count) - upper_bound);

		count += temp_count;
	}

	void miniLowPass() {
		int x_max_not_zeroed = 0;
		int y_max_not_zeroed = 0;
		
		int distance_from_centre = 0;
		for (int j = m_cols / 2; j < m_cols; ++j) {
			int i = m_rows / 2;
			double Qvalue = QM16[i][static_cast<int>(j * 16.0 / m_cols)] * 4 / (256.0 / m_data.size());

			bool zeroed = std::sqrt(std::pow(m_data[i * m_cols + j].r, 2) + std::pow(m_data[i * m_cols + j].i, 2)) / Qvalue < 1.0;
			++distance_from_centre;
			if (!zeroed) {
				x_max_not_zeroed = distance_from_centre;
			}
		}
		distance_from_centre = 0;
		for (int i = m_rows/2; i < m_rows; ++i) {
			int j = m_cols / 2;
			double Qvalue = QM16[static_cast<int>(i * 16.0 / m_rows)][j] * 4 / (256.0 / m_data.size());

			bool zeroed = std::sqrt(std::pow(m_data[i * m_cols + j].r, 2) + std::pow(m_data[i * m_cols + j].i, 2)) / Qvalue < 1.0;
			++distance_from_centre;
			if (!zeroed) {
				y_max_not_zeroed = distance_from_centre;
			}
		}
		int max_dist = std::lround(std::max(x_max_not_zeroed * m_rows / static_cast<double>(m_cols), static_cast<double>(y_max_not_zeroed))) + 3;
		double cf;
		cutOutGeneral(max_dist, 1, cf);
	}
	
	void lowPassQM() {
		int count = 0;
		int max_y_dist = 0;

		for (int i = 0; i < m_rows; ++i) {
			double Qvalue = QM16[static_cast<int>(i * 16.0 / m_rows)][static_cast<int>(i * 16.0 / m_cols)] * 4 / (256.0 / m_data.size());

			if ((std::sqrt(std::pow(m_data[i].r, 2) + std::pow(m_data[i].i, 2))) / Qvalue < 1.0) {
				max_y_dist = i - 1;
			}
		}
		for (int i = 0; i < m_cols; ++i) {
			double Qvalue = QM16[static_cast<int>(i * 16.0 / m_rows)][static_cast<int>(i * 16.0 / m_cols)] * 4 / (256.0 / m_data.size());

			if ((std::sqrt(std::pow(m_data[i].r, 2) + std::pow(m_data[i].i, 2))) / Qvalue < 1.0) {
				max_y_dist = i - 1;
			}
		}

	}

	void smartLocalisedFourier(Image& image, int startx, int starty, double& compression_ratio, double sobel, double michelson, double target_cr, bool generate_random_colour=false) {
		static int count = 0;
		if (((gradientSobelCalc(image) < sobel || michelsonContrast(image) < michelson) || image.m_data.size() <= 64) && image.m_data.size() <= 16384) {
			image.forwardFourier();
			image.varyingQM(count, target_cr);
			image.reverseFourier();
			if (!generate_random_colour) {
				for (int i = 0; i < image.m_data.size(); ++i) {
					int sectionx = i % image.m_cols;
					int sectiony = i / image.m_cols;

					int index = (starty + sectiony) * m_cols + startx + sectionx;
					m_data[index].r = image.m_data[i].r;

				}
			}
			else {
				std::random_device rd;

				std::mt19937 generator(rd());

				std::uniform_int_distribution<int> distribution(0, 255);
				int rn = distribution(generator);
				for (int i = 0; i < image.m_data.size(); ++i) {
					int sectionx = i % image.m_cols;
					int sectiony = i / image.m_cols;

					int index = (starty + sectiony) * m_cols + startx + sectionx;
					m_data[index].r = rn;
				}
			}
		}
		else {
			int width1 = image.m_cols / 2;
			int width2 = image.m_cols - width1;

			int height1 = image.m_rows / 2;
			int height2 = image.m_rows - height1;

			auto recurse = [&](int startx, int endx, int starty, int endy) {
				Image img{ endx - startx, endy - starty };



				for (int x = startx; x < endx; ++x) {
					for (int y = starty; y < endy; ++y) {
						img.m_data[(y-starty) * img.m_cols + x - startx] = image.m_data[(y) * image.m_cols + x];
					}
				}
				return img;
			};
			Image img1 = recurse(0, width2, 0, height2);
			Image img2 = recurse(width2, image.m_cols, 0, height2);
			Image img3 = recurse(0, width2, height2, image.m_rows);
			Image img4 = recurse(width2, image.m_cols, height2, image.m_rows);
			

			smartLocalisedFourier(img1, startx, starty, compression_ratio, sobel, michelson, target_cr, generate_random_colour);
			smartLocalisedFourier(img2, startx + width2, starty, compression_ratio, sobel, michelson, target_cr, generate_random_colour);
			smartLocalisedFourier(img3, startx, starty + height2, compression_ratio, sobel, michelson, target_cr, generate_random_colour);
			smartLocalisedFourier(img4, startx+width2, starty+height2, compression_ratio, sobel, michelson, target_cr, generate_random_colour);
		}

		if (image.m_data.size() == m_data.size()) {
			compression_ratio = count / static_cast<double>(m_data.size());
			count = 0;
		}


	}

	void FFTShift() {
		int half_row = m_rows / 2;
		int half_col = m_cols / 2;
		std::vector<kiss_fft_cpx> temp_data = m_data;

		for (int i = 0; i < m_rows; ++i) {
			for (int j = 0; j < m_cols; ++j) {
				if (j >= half_col && i >= half_row) {
					m_data[i * m_cols + j] = temp_data[(m_rows + half_row - i-1) * m_cols + (m_cols - j + half_col - 1)];
				}

				if (j < half_col && i >= half_row) {
					m_data[i * m_cols + j] = temp_data[(m_rows + half_row - i-1) * m_cols + (half_col - j - 1)];
				}

				if (j >= half_col && i < half_row) {
					m_data[i * m_cols + j] = temp_data[(half_row - i - 1) * m_cols + m_cols - j + half_col-1];
				}

				if (j < half_col && i < half_row) {
					m_data[i* m_cols + j] = temp_data[(half_row - i - 1) * m_cols + half_col -j - 1];
				}
			}
		}
	}

	double getMean() const {
		double mean = 0;
		for (int i = 0; i < m_data.size(); ++i) {
			mean += m_data[i].r;
		}
		mean /= m_data.size();
		return mean;
	}

	double getSD(double mean) const {
		double variance = 0;
		double standard_deviation = 0;

		for (int i = 0; i < m_data.size(); ++i) {
			variance += std::pow((m_data[i].r - mean), 2);
		}
		variance /= m_data.size();
		standard_deviation = std::sqrt(variance);
		return standard_deviation;
	}

	int getRows() const { return m_rows; }
	int getCols() const { return m_cols; }
	int getPixelCount() const { return m_data.size(); }
	std::string getName() const { return m_name; }

	void setData(const Image& image) {
		m_data = image.m_data;
	}

	friend double calculateSSIM(Image& image1, Image& image2);

private:
	std::string m_path;
	std::string m_name;

	int m_rows, m_cols;

	std::vector<kiss_fft_cpx> m_data;

	const std::array<std::array<int, 8>, 8> default_QM{ {
		{16, 11, 10, 16, 24, 40, 51, 61},
		{12, 12, 14, 19, 26, 58, 60, 55},
		{14, 13, 16, 24, 40, 57, 69, 56},
		{14, 17, 22, 29, 51, 87, 80, 62},
		{18, 22, 37, 56, 68, 109, 103, 77},
		{24, 35, 55, 64, 81, 104, 113, 92},
		{49, 64, 78, 87, 103, 121, 120, 101},
		{72, 92, 95, 98, 112, 100, 103, 99}
	} };

	const std::array<std::array<int, 16>, 16> QM16 { {
		{99, 103, 100, 112, 98, 95, 92, 72, 72, 92, 95, 98, 112, 100, 103, 99},
		{101, 120, 121, 103, 87, 80, 62, 51, 51, 62, 80, 87, 103, 121, 120, 101},
		{92, 113, 104, 81, 64, 55, 35, 24,24, 35, 55, 64, 81, 104, 113, 92 },
		{77, 103, 109, 68, 56, 37, 22, 18,18, 22, 37, 56, 68, 109, 103, 77 },
		{62, 80, 87, 51, 29, 22, 17, 14,14, 17, 22, 29, 51, 87, 80, 62 },
		{56, 69, 57, 40, 24, 16, 13, 14,14, 13, 16, 24, 40, 57, 69, 56 },
		{55, 60, 58, 26, 19, 14, 12, 12,12, 12, 14, 19, 26, 58, 60, 55 },
		{61, 51, 40, 24, 16, 10, 11, 16,16, 11, 10, 16, 24, 40, 51, 61 },
		{ 61, 51, 40, 24, 16, 10, 11, 16, 16, 11, 10, 16, 24, 40, 51, 61 },
		{55, 60, 58, 26, 19, 14, 12, 12, 12, 12, 14, 19, 26, 58, 60, 55},
		{ 56, 69, 57, 40, 24, 16, 13, 14, 14, 13, 16, 24, 40, 57, 69, 56 },
		{ 62, 80, 87, 51, 29, 22, 17, 14,14, 17, 22, 29, 51, 87, 80, 62  },
		{ 77, 103, 109, 68, 56, 37, 22, 18,18, 22, 37, 56, 68, 109, 103, 77 },
		{ 92, 113, 104, 81, 64, 55, 35, 24,24, 35, 55, 64, 81, 104, 113, 92 },
		{ 101, 120, 121, 103, 87, 80, 62, 51,51, 62, 80, 87, 103, 121, 120, 101 },
		{ 99, 103, 100, 112, 98, 95, 92, 72,72, 92, 95, 98, 112, 100, 103, 99 }
	}};
};

inline double calculateSSIM(Image& image1, Image& image2) {
	if (image1.m_data.size() != image2.m_data.size()) {
		std::cerr << "Image Size Mismatch when Calculating SSIM!\n";
		return 1;
	}

	// Calculate Covariance
	double covariance = 0;
	double L = 255;
	int image_size = image1.m_data.size();

	double mean1 = image1.getMean();
	double mean2 = image2.getMean();

	double SD1 = image1.getSD(mean1);
	double SD2 = image2.getSD(mean2);

	for (int i = 0; i < image_size; ++i) {
		covariance += (image1.m_data[i].r - mean1) * (image2.m_data[i].r - mean2);
	}
	covariance /= image_size;

	double C1 = std::pow((0.01 * L), 2);
	double C2 = std::pow((0.03 * L), 2);
	
	double luminance = (2 * mean1 * mean2 + C1) / (mean1 * mean1 + mean2 * mean2 + C1);
	double structure = (covariance + C2 / 2) / (SD1 * SD2 + C2 / 2);
	double contrast = (2 * SD1 * SD2 + C2) / (SD1 * SD1 + SD2 * SD2 + C2);

	return luminance * structure * contrast;
}
