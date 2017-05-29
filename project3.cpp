// Capturing Reality Projekt 2
// HDR Recovery von Projekt 1
// Blending mit Laplace Pyramide (mithilfe von OpenCV für leichteres filtern)
// von Janina Hüther, Lennart Jarms
#include "CImg.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
//#include "opencv2/opencv.hpp"

# define M_PI 3.14159265358979323846

using namespace std;
using namespace cimg_library;

typedef struct ImageInfo {
	std::string imageName;
	float t;
};

typedef struct rgb {
	float r;
	float g;
	float b;

	rgb() :r(0.f),
		g(0.f),
		b(0.f) {}
	rgb(float _r, float _g, float _b) :r(_r),
		g(_g),
		b(_b) {}
};

// =============
// = Project 1 =
// =============

vector<ImageInfo> readFile(string filename) {
	ifstream in(filename);
	string line;
	if (!in.is_open()) {
		cout << "loadHDRGEN: can not open " << endl;
	}

	vector<ImageInfo> imgVector;
	string img;
	float t;
	const int MAX = 256;
	for (string img; in >> setw(MAX) >> img; ) {
		in >> setw(MAX) >> t;
		ImageInfo info;
		info.imageName = img;
		info.t = t;
		imgVector.push_back(info);
	}

	in.close();

	return imgVector;
}

//reading csv file http://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c
class CSVRow
{
public:
	std::string const& operator[](std::size_t index) const
	{
		return m_data[index];
	}
	std::size_t size() const
	{
		return m_data.size();
	}
	void readNextRow(std::istream& str)
	{
		std::string         line;
		std::getline(str, line);

		std::stringstream   lineStream(line);
		std::string         cell;

		m_data.clear();
		while (std::getline(lineStream, cell, ','))
		{
			m_data.push_back(cell);
		}
		// This checks for a trailing comma with no data after it.
		if (!lineStream && cell.empty())
		{
			// If there was a trailing comma then add an empty element.
			m_data.push_back("");
		}
	}
private:
	std::vector<std::string>    m_data;
};

std::istream& operator >> (std::istream& str, CSVRow& data)
{
	data.readNextRow(str);
	return str;
}

vector<rgb> readResponseCurve() {
	vector<rgb> I;
	std::ifstream file("responseCurve_precalculated.csv");

	CSVRow row;
	while (file >> row)
	{
		I.push_back({ std::stof(row[1]),std::stof(row[2]),std::stof(row[3]) });
	}

	return I;
}

//schreibt response curve in die Datei "responseCurve.csv" im build ordner
void writeToCSV(vector<rgb> I) {
	ofstream myfile;
	string r = "", g = "", b = "";
	myfile.open("responseCurve.csv");
	for (int i = 0; i < I.size(); i++) {
		r = std::to_string(I[i].r);
		g = std::to_string(I[i].g);
		b = std::to_string(I[i].b);
		myfile << std::to_string(i) << "," << r << "," << g << "," << b << "\n";
	}
	myfile.close();
	return;
}

//Hilfsfunktion zum Aufrufen von Gnuplot
void plotGnuPlot() {
	system("gnuplot -p -e \"load 'responseCurve.p'\"");
}

//w
float weight(float y) {
	if (y <= 13.f || y >= 250.f) return 0;
	else return std::exp(-4.f * (powf(y - 127.5f, 2) / powf(127.5f, 2)));
}

//xj
CImg<float> calculate_irradiance(vector<CImg<float>> images, vector<ImageInfo> imgInfo, vector<rgb> I) {
	CImg<float> result(images[0]); //bild mit gleichen dimensionen anlegen
	result.fill(0);
	vector<float> numerator(images[0].size(), 0.0f);
	vector<float> sum(images[0].size(), 0.0f);
	float Iy, t, w;
	CImg<float> img;

	for (int i = 0; i < images.size(); i++) {
		img = images[i];
		t = 1.0f / imgInfo[i].t; //exposure time

		cimg_foroff(img, j) { //iteriert über bild buffer
			w = weight(img[j]);

			//response function value for pixel value
			if (j < img.size() / 3) { //Red values
				Iy = I[img[j]].r;
			}
			else if (j < 2 * img.size() / 3) { //Green values 
				Iy = I[img[j]].g;
			}
			else { //Blue values
				Iy = I[img[j]].b;
			}

			numerator[j] += (w * t * Iy);
			sum[j] += (w * powf(t, 2));
		}
	}

	for (int y = 0; y < result.size(); y++) {
		if (sum[y] != 0.0f) {
			result[y] = numerator[y] / sum[y];
		}
		else {
			result[y] = 255;
			//result[y] = 0;
		}
		result[y] = result[y] > 255 ? 255 : result[y] < 0 ? 0 : result[y];
	}

	return result;
}

//I
vector<rgb> calculate_response_curve(vector<CImg<float>> images, vector<ImageInfo> imgInfo, CImg<float> x) {
	vector<rgb> I(256, { 0.f,0.f,0.f });
	CImg<float> img;
	float t;

	for (int m = 0; m < I.size(); m++) {
		rgb Card = { 0.f,0.f,0.f };
		for (int i = 0; i < images.size(); i++) {
			img = images[i];
			t = 1.0f / imgInfo[i].t;
			cimg_foroff(img, j) { //iteriert über bild buffer
				if ((int)img[j] == m) {
					if (j < img.size() / 3) { //Red values
						I[m].r += t * x[j];
						Card.r += 1;
					}
					else if (j < 2 * img.size() / 3) { //Green values 
						I[m].g += t * x[j];
						Card.g += 1;
					}
					else { //Blue values
						I[m].b += t * x[j];
						Card.b += 1;
					}
				}
			}
		}
		I[m].r /= Card.r;
		I[m].g /= Card.g;
		I[m].b /= Card.b;
	}

	//normalisieren sodass I_128 = 1.0f ist
	rgb i128 = I[128];
	for (int i = 0; i < 256; i++) {
		I[i].r /= i128.r;
		I[i].g /= i128.g;
		I[i].b /= i128.b;
	}

	return I;
}

//O
rgb calculate_objective_f(vector<CImg<float>> images, vector<ImageInfo> imgInfo, vector<rgb> I, CImg<float> x) {
	rgb result = { 0.f,0.f,0.f };
	float t;
	CImg<float> img;

	for (int i = 0; i < images.size(); i++) {
		t = 1.0f / imgInfo[i].t;
		img = images[i];
		cimg_foroff(img, j) { //iteriert über bild buffer
			if (j < img.size() / 3) { //Red values
				result.r += weight(img[j]) * (powf(I[img[j]].r - (t*x[j]), 2));
			}
			else if (j < 2 * img.size() / 3) { //Green values 
				result.g += weight(img[j]) * (powf(I[img[j]].g - (t*x[j]), 2));
			}
			else { //Blue values
				result.b += weight(img[j]) * (powf(I[img[j]].g - (t*x[j]), 2));
			}
		}
	}

	return result;
}

//Tonemapping von Drago et al. 2003
CImg<float> calculate_tone_mapping(CImg<float> x) {
	CImg<float> result_xyz = x.get_RGBtoXYZ();
	float b = 0.75f; //0.0 - 1.0
	float L_avg = 0.f;
	float Lw_max = 0.f;
	float Ld_max = 100.f;

	cimg_forXY(result_xyz, x, y) {
		L_avg += result_xyz(x, y, 0, 1);
		if (result_xyz(x, y, 0, 1) > Lw_max)
			Lw_max = result_xyz(x, y, 0, 1);
	}
	int w = result_xyz.width(), h = result_xyz.height();
	L_avg /= (w * h);

	Lw_max /= L_avg; //max luminanz normalisieren

	CImg<float> tmpImg(result_xyz);

	cimg_forXY(result_xyz, x, y) {
		float L_w = result_xyz(x, y, 0, 1) / L_avg;
		float L_d_left = (Ld_max * 0.01f) / log10(Lw_max + 1.f);
		float L_d_right = (log(L_w + 1.f) / log(2.f + powf(L_w / Lw_max, log(b) / log(.5f)) * 8.f));
		float L_d = L_d_left * L_d_right;

		tmpImg(x, y, 0, 1) = L_d;
	}

	//skalieren der channel 0 und 2 abh. von channel 1
	cimg_forXY(result_xyz, x, y) {
		float s = tmpImg(x, y, 0, 1) / result_xyz(x, y, 0, 1);
		result_xyz(x, y, 0, 0) *= s;
		result_xyz(x, y, 0, 1) *= s;
		result_xyz(x, y, 0, 2) *= s;
	}

	return result_xyz.get_XYZtoRGB();
}

CImg<float> createHDRImage(string name, bool loadResponseCurve) {
	//.hdrgen einlesen
	vector<ImageInfo> imageNames = readFile("img/HDRsequence/" + name);
	for (std::vector<ImageInfo>::const_iterator i = imageNames.begin(); i != imageNames.end(); ++i)
		std::cout << (*i).imageName << " t: " << 1. / (*i).t << "\n";


	//Alle bilder einlesen
	std::vector<CImg<float>> images;
	for (int i = 0; i < imageNames.size(); i++) {
		std::string imageName = "img/HDRsequence/" + imageNames[i].imageName;
		images.push_back(CImg<>(imageName.c_str()));
	}

	CImg<float> x;
	if (loadResponseCurve == true) {
		vector<rgb> I = readResponseCurve();
		x = calculate_irradiance(images, imageNames, I);
	}
	else {
		//initiale response curve mit I_1 = 1/128, I_128 = 1, I_256 = 2f, 
		vector<rgb> I(256);
		for (int i = 0; i < 256; i++) {
			float v = (i) / 128.f;
			I[i] = { v,v,v };
		}
		x = calculate_irradiance(images, imageNames, I);
		cout << "Initiale Iteration: " << endl;
		rgb o = calculate_objective_f(images, imageNames, I, x);
		cout << "O.rgb : " << o.r << ", " << o.g << ", " << o.b << endl;

		//TODO convergenz durch delta oder so in while schleife abfragen
		rgb o_new;
		bool converged;
		int n = 1;
		do {
			cout << "Iteration: " << n << endl;
			I = calculate_response_curve(images, imageNames, x);
			x = calculate_irradiance(images, imageNames, I);
			o_new = calculate_objective_f(images, imageNames, I, x);

			converged =
				(o.r - o_new.r) / o.r < 0.1f &&
				(o.g - o_new.g) / o.g < 0.1f &&
				(o.b - o_new.b) / o.b < 0.1f;

			n += 1;
			o = o_new;
			cout << "O.rgb : " << o.r << ", " << o.g << ", " << o.b << endl;
		} while (!converged);
		writeToCSV(I);
	}

	
	CImg<float> image = calculate_tone_mapping(x);

	return image;
}

// =============
// = Project 2 =
// =============

float
BilinearInterpolation(float q11, float q12, float q21, float q22, float x1, float x2, float y1, float y2, float x, float y)
{
	float R1, R2, P;
	if (x1 == x2) {
		R1 = (x2 - x) *q11;
		R2 = (x2 - x)*q12;
	}
	else {
		R1 = ((x2 - x) / (x2 - x1))*q11 + ((x - x1) / (x2 - x1))*q21;
		R2 = ((x2 - x) / (x2 - x1))*q12 + ((x - x1) / (x2 - x1))*q22;
	}

	if (y1 == y2) {
		P = (y2 - y)*R1;
	}
	else {
		P = ((y2 - y) / (y2 - y1))*R1 + ((y - y1) / (y2 - y1))*R2;
	}

	P = P == 0 ? q12 : P;

	return P;
}

//Mit bilinearer Interpolation
CImg<float> backwardWarp(CImg<float> src) {
	CImg<> image(src.width()*2, src.height()*0.5f, 1, 3);
	image.fill(0.f);
	float _x = 0.f, _y = 0.f, _z = 0.f;
	
	float radius = src.width() < src.height() ? src.width() / 2.f : src.height() / 2.f;

	cimg_forXY(image, x, y) {
		//Screen to cartesian
		_y = (src.height() - y);

		//normalize
		_x = (x / radius) - 1.f;
		_y = (_y / radius) - 1.f;

		////Cart to cyl
		//float phi_cyl = M_PI * _x;
		float theta_cyl = M_PI * _y;

		////Cyl to sphere
		float theta_sphere = theta_cyl;
		float phi_sphere = atan(_x* M_PI) ;
		//Sphere to image coords

		float phi = phi_sphere;		
		float theta = theta_sphere;

		float phi1 = atan2(cos(phi) * sin(theta), cos(theta));
		float theta1 = atan2(sqrt(1.0 - powf(sin(phi), 2) * powf(sin(theta), 2)), sin(phi) * sin(theta));
		phi = phi1;
		theta = theta1;
		phi1 = atan2(cos(phi) * sin(theta), cos(theta));
		theta1 = atan2(sqrt(1.0 - powf(sin(phi), 2) * powf(sin(theta), 2)), sin(phi) * sin(theta));
		/*float phi1 = atan2(cos(theta), sin(phi)*sin(theta));
		float theta1 = atan2(sqrt(1.0 + powf(sin(theta), 2) * powf(sin(phi), 2))- powf(sin(theta), 2), (cos(phi)*sin(theta)) );*/
		float r = sin(theta1 / 2.f);
		float xi = r * cos(phi1);
		float yi = r * sin(phi1);


		//for bilinear interpolation test (stretches original img)
		//_x = _x / 2.f;
		//float xi = (_x + 0.5f)*radius;
		//float yi = (_y + 1)*radius;

		xi = (xi + 1)*radius;
		yi = (yi + 1)*radius;
		//bilinear interpolation with replicate boundary conditions
		float x1 = floor(xi) >= src.width() ? src.width() - 1 : floor(xi);
		float x2 = ceil(xi) >= src.width() ? src.width() - 1 : ceil(xi);
		float y2 = floor(yi) >= src.height() ? src.height() - 1 : floor(yi);
		float y1 = ceil(yi) >= src.height() ? src.height() - 1 : ceil(yi);
		//cout << "_x, _y: " << _x << " " << _y << endl;
		//cout << "x, y: " << x << " " << y << endl;

		//cout << "x1, x2: " << x1 << " " <<  x2 << " y1, y2: " << y1  << " " << y2 << endl;


		float q12_r = src(x1, y1, 0), q12_g = src(x1, y1, 1), q12_b = src(x1, y1, 2);
		float q22_r = src(x1, y2, 0), q22_g = src(x1, y2, 1), q22_b = src(x1, y2, 2);
		float q11_r = src(x2, y1, 0), q11_g = src(x2, y1, 1), q11_b = src(x2, y1, 2);
		float q21_r = src(x2, y2, 0), q21_g = src(x2, y2, 1), q21_b = src(x2, y2, 2);

		//cout << "q11: " << q11_r << ", " << "q12: " << q12_r << ", " << "q21: " << q21_r << ", " << "q22: " << q22_r << ", " << " x1, x2, y1, y2: " << x1 << ", " << x2 << ", " << y1 << ", " << y2 << ", " << " xi, yi:" << xi << ", " << yi  << endl;

		float _r = BilinearInterpolation(q11_r, q12_r, q21_r, q22_r, x1, x2, y1, y2, xi, yi);
		float _g = BilinearInterpolation(q11_g, q12_g, q21_g, q22_g, x1, x2, y1, y2, xi, yi);
		float _b = BilinearInterpolation(q11_b, q12_b, q21_b, q22_b, x1, x2, y1, y2, xi, yi);
		//undo normalize

		//cout << "r,g,b: " << _r << ", " << _g << ", " << _b << endl;
		//y = src.height() - y;
		image(x, y, 0) = _r;
		image(x, y, 1) = _g;
		image(x, y, 2) = _b;
	}

	return image;
}

vector<CImg<float>> createGaussianPyramid(const CImg<float> image, int n_level) {
	vector<CImg<float>> gpyramid(n_level);
	gpyramid[0] = image;
	int width, height;
	for (int i = 1; i < n_level; i++) {
		width = gpyramid[i - 1].width() / 2;
		height = gpyramid[i - 1].height() / 2;
		gpyramid[i] = gpyramid[i - 1].get_blur(3.f).get_resize(width, height, 1, image.spectrum(), 3);
	}

	return gpyramid;
}

vector<CImg<float>> createLaplacianPyramid(const CImg<float> image, int n_level) {
	vector<CImg<float>> lpyramid(n_level);
	vector<CImg<float>> gpyramid = createGaussianPyramid(image, n_level);
	int width, height;
	for (int i = 0; i < n_level - 1; i++) {
		width = gpyramid[i].width();
		height = gpyramid[i].height();
		lpyramid[i] = gpyramid[i];

		CImg<float> gpyr = gpyramid[i + 1].get_resize(width, height, 1, 3, 3);
		cimg_forXYC(lpyramid[i], x, y, c) {
			float val = (gpyramid[i](x, y, c) - gpyr(x, y, c));
			val = val > 255 ? 255 : val < 0 ? 0 : val;
			lpyramid[i](x, y, c) = val;
		}
		//lpyramid[i] = gpyramid[i] - gpyramid[i + 1].get_resize(width, height, 1, 3, 3);

		if (i == 0) lpyramid[i] *= 1.5f;
		else if (i == 1) lpyramid[i] *= 1.9f;
	}
	lpyramid[n_level - 1] = gpyramid[n_level - 1];
	return lpyramid;
}

//Binärmaske zum Blenden
vector<CImg<float>> createBlendingMask(int width, int left_width, int height, int n_level, float overlap) {
	CImg<float> mask = CImg<float>(width, height, 1, 1, 0);
	cimg_forXY(mask, x, y) {
		if (x  < overlap * left_width)
			mask(x, y) = 1;
	}
	vector<CImg<float>> gpyramid_mask = createGaussianPyramid(mask, n_level);
	return gpyramid_mask;
}

CImg<float> reconstructImage(vector<CImg<float>> lpyramid, int n_level) {
	CImg<float> image = lpyramid[n_level - 1];
	int width, height;
	for (int i = n_level - 2; i >= 0; i--) {
		width = lpyramid[i].width();
		height = lpyramid[i].height();
		image.resize(width, height, 1, 3, 3);

		cimg_forXYC(lpyramid[i], x, y, c) {
			float val = (image(x, y, c) + lpyramid[i](x, y, c));
			val = val > 255 ? 255 : val < 0 ? 0 : val;
			image(x, y, c) = val;
		}
	}

	return image;
}

//gibt geblendete laplace pyramiden zurück
CImg<float> blendImages(const CImg<float> left, const CImg<float> right) {
	int min_length = (left.width() < left.height()) ? left.width() : left.height();
	int n_level = floor(log2(min_length)) - 5; //kleinst mögliches pyramiden level

	vector<CImg<float>> l_lpyramid = createLaplacianPyramid(left, n_level); //Laplace Pyramide linkes bild
	vector<CImg<float>> r_lpyramid = createLaplacianPyramid(right, n_level); //Laplace Pyramide rechtes bild
	
	//Wir wählen eine feste Nahstelle bei 75% des Bildes
	//0.42 und 0.45 für entfernen der kamera
	// 0.55 und 0.55 für doppel kamera
	// 0.38 und 0.75 für p1, p2
	float overlap_left = 0.65;
	float overlap_right = 0.65;
	//int new_width = 2 * overlap * left.width(); //
	int new_width = overlap_left * left.width() + overlap_right * right.width(); //
	cout << "new width: " << new_width << endl;
	vector<CImg<float>> mask = createBlendingMask(new_width, left.width(), left.height(), n_level, overlap_left);

	vector<CImg<float>> blended_lpyramid(n_level);
	for (int i = 0; i < n_level; i++) {
		new_width = overlap_left*l_lpyramid[i].width() + overlap_right * r_lpyramid[i].width();
		blended_lpyramid[i] = CImg<float>(new_width, l_lpyramid[i].height(), 1, l_lpyramid[i].spectrum(), 0);
		//blended_lpyramid[i] = CImg<float>(l_lpyramid[i].width(), l_lpyramid[i].height(), 1, 3);
		cout << "mask w,h: " << mask[i].width() << ", " << mask[i].height() << endl;
		cimg_forXYC(blended_lpyramid[i], x, y,c) {
			//TODO abfrage
			if (x < overlap_left * (l_lpyramid[i].width())) {
				if (mask[i](x, y) < 1.f) {
					blended_lpyramid[i](x, y, c) = mask[i](x, y) * l_lpyramid[i](x, y, c) + (1- mask[i](x, y)) * r_lpyramid[i](x - overlap_left * l_lpyramid[i].width() + (1 - overlap_right) * r_lpyramid[i].width(), y, c);
				}
				else {
					blended_lpyramid[i](x, y, c) = l_lpyramid[i](x, y, c);
				}
				
			}
			else {
				if (mask[i](x, y) > 0.f) {
					blended_lpyramid[i](x, y, c) = mask[i](x, y) * l_lpyramid[i](x, y, c) + (1 - mask[i](x, y)) * r_lpyramid[i](x - overlap_left * l_lpyramid[i].width() + (1 - overlap_right) * r_lpyramid[i].width(), y, c);
				}
				else {
					blended_lpyramid[i](x, y, c) = r_lpyramid[i](x - overlap_left * l_lpyramid[i].width() + (1 - overlap_right) * r_lpyramid[i].width(), y, c);
				}
			}
			
			//blended_lpyramid[i](x, y, c) =mask[i](x, y) == 1 ? l_lpyramid[i](x, y, c) : r_lpyramid[i](x, y, c);
		}
	}

	CImg<float> blended_image = reconstructImage(blended_lpyramid, n_level);

	return blended_image;
}

int main(int argc, char **argv) {
	/*string imageName = "img/spherical_panorama.ppm";
	CImg<> src(imageName.c_str());*/
	CImg<> image;
	//string imageName = "img/HDRsequence/a0s.ppm";
	string imageName1 = "img/p0.ppm";
	string imageName2 = "img/p2.ppm";
	string hdrFile = "e.hdrgen";
	bool loadResponseCurve = true; //berechnete response curve verwenden
	//CImg<> src = createHDRImage(hdrFile, loadResponseCurve);
	CImg<> src1(imageName1.c_str());
	CImg<> src2(imageName2.c_str());

	CImg<float> image_warped1 = backwardWarp(src1); //ohne
	CImg<float> image_warped2 = backwardWarp(src2); //ohne

	//CImg<float> blended_image = blendImages(image_warped1, image_warped2);
	//image = blended_image;

	image = image_warped1;
	
	// ==============================
	// = Code vom CImg tutorial.cpp =
	// ==============================
	/*
	#
	#  File        : tutorial.cpp
	#                ( C++ source file )
	#
	#  Description : View the color profile of an image, along the X-axis.
	#                This file is a part of the CImg Library project.
	#                ( http://cimg.eu )
	#
	#  Copyright   : David Tschumperle
	#                ( http://tschumperle.users.greyc.fr/ )
	#
	#  License     : CeCILL v2.0
	#                ( http://www.cecill.info/licences/Licence_CeCILL_V2-en.html )
	#
	#  This software is governed by the CeCILL  license under French law and
	#  abiding by the rules of distribution of free software.  You can  use,
	#  modify and/ or redistribute the software under the terms of the CeCILL
	#  license as circulated by CEA, CNRS and INRIA at the following URL
	#  "http://www.cecill.info".
	#
	#  As a counterpart to the access to the source code and  rights to copy,
	#  modify and redistribute granted by the license, users are provided only
	#  with a limited warranty  and the software's author,  the holder of the
	#  economic rights,  and the successive licensors  have only  limited
	#  liability.
	#
	#  In this respect, the user's attention is drawn to the risks associated
	#  with loading,  using,  modifying and/or developing or reproducing the
	#  software by the user in light of its specific status of free software,
	#  that may mean  that it is complicated to manipulate,  and  that  also
	#  therefore means  that it is reserved for developers  and  experienced
	#  professionals having in-depth computer knowledge. Users are therefore
	#  encouraged to load and test the software's suitability as regards their
	#  requirements in conditions enabling the security of their systems and/or
	#  data to be ensured and,  more generally, to use and operate it in the
	#  same conditions as regards security.
	#
	#  The fact that you are presently reading this means that you have had
	#  knowledge of the CeCILL license and that you accept its terms.
	#
	*/
	// Create two display window, one for the image, the other for the color profile.
	CImgDisplay
		main_disp(image, "Color image", 0),
		draw_disp(500, 400, "Color profile of the X-axis", 0);

	// Define colors used to plot the profile, and a hatch to draw the vertical line
	unsigned int hatch = 0xF0F0F0F0;
	const unsigned char
		red[] = { 255,0,0 },
		green[] = { 0,255,0 },
		blue[] = { 0,0,255 },
		black[] = { 0,0,0 };



	// Enter event loop. This loop ends when one of the two display window is closed or
	// when the keys 'ESC' or 'Q' are pressed.
	while (!main_disp.is_closed() && !draw_disp.is_closed() &&
		!main_disp.is_keyESC() && !draw_disp.is_keyESC() && !main_disp.is_keyQ() && !draw_disp.is_keyQ()) {

		// Handle display window resizing (if any)
		if (main_disp.is_resized()) main_disp.resize().display(image);
		draw_disp.resize();

		if (main_disp.mouse_x() >= 0 && main_disp.mouse_y() >= 0) { // Mouse pointer is over the image

			const int
				xm = main_disp.mouse_x(),                     // X-coordinate of the mouse pointer over the image
				ym = main_disp.mouse_y(),                     // Y-coordinate of the mouse pointer over the image
				xl = xm*draw_disp.width() / main_disp.width(),  // Corresponding X-coordinate of the hatched line
				x = xm*image.width() / main_disp.width(),     // Corresponding X-coordinate of the pointed pixel in the image
				y = ym*image.height() / main_disp.height();   // Corresponding Y-coordinate of the pointex pixel in the image

															  // Retrieve color component values at pixel (x,y)
			const unsigned int
				val_red = image(x, y, 0),
				val_green = image(x, y, 1),
				val_blue = image(x, y, 2);

			// Create and display the image of the intensity profile
			CImg<unsigned char>(draw_disp.width(), draw_disp.height(), 1, 3, 255).
				draw_grid(-50 * 100.0f / image.width(), -50 * 100.0f / 256, 0, 0, false, true, black, 0.2f, 0xCCCCCCCC, 0xCCCCCCCC).
				draw_axes(0, image.width() - 1.0f, 255.0f, 0.0f, black).
				draw_graph(image.get_shared_row(y, 0, 0), red, 1, 1, 0, 255, 1).
				draw_graph(image.get_shared_row(y, 0, 1), green, 1, 1, 0, 255, 1).
				draw_graph(image.get_shared_row(y, 0, 2), blue, 1, 1, 0, 255, 1).
				draw_text(30, 5, "Pixel (%d,%d)={%d %d %d}", black, 0, 1, 16,
					main_disp.mouse_x(), main_disp.mouse_y(), val_red, val_green, val_blue).
				draw_line(xl, 0, xl, draw_disp.height() - 1, black, 0.5f, hatch = cimg::rol(hatch)).
				display(draw_disp);
		}
		else
			// else display a text in the profile display window.
			CImg<unsigned char>(draw_disp.width(), draw_disp.height()).fill(255).
			draw_text(draw_disp.width() / 2 - 130, draw_disp.height() / 2 - 5, "Mouse pointer is outside the image",
				black, 0, 1, 16).display(draw_disp);

		// Temporize event loop
		cimg::wait(20);
	}

	return 0;
}