#pragma once

#include <cstddef>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <type_traits>
#include <iostream>

struct imgdta {
	imgdta(std::size_t width, std::size_t height = 1, std::size_t channels = 1, std::size_t bps = 1, bool is_float = false) : w(width), h(height), nc(channels), nb(bps), _data(w * h * nc * bps), flt(is_float)
	{}
	imgdta() : imgdta(0, 0)
	{}

	std::size_t width() const
	{
		return w;
	}
	std::size_t height() const
	{
		return h;
	}
	std::size_t channels() const
	{
		return nc;
	}
	std::size_t bps() const
	{
		return nb;
	}
	std::size_t pixels() const
	{
		return width() * height();
	}
	std::size_t size() const
	{
		return pixels() * channels();
	}
	std::size_t bytes() const
	{
		return size() * bps();
	}
	std::size_t stride() const
	{
		return width() * channels();
	}
	std::size_t pitch() const
	{
		return stride() * bps();
	}
	bool is_float() const
	{
		return flt;
	}

	void resize(std::size_t width, std::size_t height = 1, std::size_t channels = 1, std::size_t bps = 1)
	{
		w = width;
		h = height;
		nc = channels;
		nb = bps;
		_data.resize(w * h * nc * nb);
	}

	// raw indexing, ignoring channels
	uint8_t *data(std::size_t idx = 0)
	{
		return _data.data() + idx * bps();
	}
	const uint8_t *data(std::size_t idx = 0) const
	{
		return _data.data() + idx * bps();
	}
	// 1D indexing of the image matrix
	uint8_t *ptr(std::size_t idx, std::size_t c = 0)
	{
		return data(idx * nc + c);
	}
	const uint8_t *ptr(std::size_t idx, std::size_t c = 0) const
	{
		return data(idx * nc + c);
	}

	// 2D indexing of the image matrix
	uint8_t *ptr2d(std::size_t x, std::size_t y, std::size_t c = 0)
	{
		return data(y * w * nc + x * nc + c);
	}
	const uint8_t *ptr2d(std::size_t x, std::size_t y, std::size_t c = 0) const
	{
		return data(y * w * nc + x * nc + c);
	}

private:
	std::size_t w, h, nc, nb;
	bool flt;
	std::vector<uint8_t> _data;
};

template <typename T>
struct image : imgdta {
	image(const imgdta &img) : imgdta(img)
	{}
	image(std::size_t width, std::size_t height = 1, std::size_t channels = 1) : imgdta(width, height, channels, sizeof(T), std::is_floating_point<T>::value)
	{}
	image() : image(0, 0)
	{}
	void resize(std::size_t width, std::size_t height = 1, std::size_t channels = 1)
	{
		imgdta::resize(width, height, channels, sizeof(T));
	}

	// raw indexing, ignoring channels
	T *data(std::size_t idx = 0)
	{
		return (T*)imgdta::data(idx);
	}
	const T *data(std::size_t idx = 0) const
	{
		return (const T*)imgdta::data(idx);
	}

	// 1D indexing of the image matrix
	T *ptr(std::size_t idx, std::size_t c = 0)
	{
		return data(idx * channels() + c);
	}
	const T *ptr(std::size_t idx, std::size_t c = 0) const
	{
		return data(idx * channels() + c);
	}

	// 2D indexing of the image matrix
	T *ptr2d(std::size_t x, std::size_t y, std::size_t c = 0)
	{
		return data(y * width() * channels() + x * channels() + c);
	}
	const T *ptr2d(std::size_t x, std::size_t y, std::size_t c = 0) const
	{
		return data(y * width() * channels() + x * channels() + c);
	}

	// raw indexing of the image matrix
	T &operator[](std::size_t idx)
	{
		return *data(idx);
	}
	const T &operator[](std::size_t idx) const
	{
		return *data(idx);
	}

	// 1D indexing of the image matrix
	T &at(std::size_t idx, std::size_t c = 0)
	{
		return *ptr(idx, c);
	}
	const T &at(std::size_t idx, std::size_t c = 0) const
	{
		return *ptr(idx, c);
	}

	// 2D indexing of the image matrix
	T &at2d(std::size_t x, std::size_t y, std::size_t c = 0)
	{
		return *ptr2d(x, y, c);
	}
	const T &at2d(std::size_t x, std::size_t y, std::size_t c = 0) const
	{
		return *ptr2d(x, y, c);
	}
	T &operator()(std::size_t x, std::size_t y, std::size_t c = 0)
	{
		return *ptr2d(x, y, c);
	}
	const T &operator()(std::size_t x, std::size_t y, std::size_t c = 0) const
	{
		return *ptr2d(x, y, c);
	}

	template <typename U = float>
	U at_lin(float x, float y, std::size_t c = 0) const
	{
		std::size_t xl = x, yl = y;
		std::size_t xu = xl + 1, yu = yl + 1;
		float ax = x - xl, ay = y - yl;

		// the user has to take care about the other boundaries
		if (xu == width()) xu = width() - 1;
		if (yu == height()) yu = height() - 1;

		U p0 = (*this)(xl, yl, c), p1 = (*this)(xu, yl, c), p2 = (*this)(xl, yu, c), p3 = (*this)(xu, yu, c);
		U l0 = p0 * (1.f - ax) + p1 * ax;
		U l1 = p2 * (1.f - ax) + p3 * ax;
		return l0 * (1.f - ay) + l1 * ay;
	}

	// vector variants
	template <typename V>
	T &at_v(const V &v, std::size_t c = 0)
	{
		return (*this)(v[0], v[1], c);
	}
	template <typename V>
	const T &at_v(const V &v, std::size_t c = 0) const
	{
		return (*this)(v[0], v[1], c);
	}
	template <typename U = float, typename V>
	U at_v_lin(const V &v, std::size_t c = 0) const
	{
		return at_lin(v[0], v[1], c);
	}
};

typedef image<uint8_t> image_b;
typedef image<uint16_t> image_s;
typedef image<float> image_f;

namespace image_io {

#ifdef WITH_PNG
imgdta load_png(const char *filename);
void save_png(const imgdta &image, const char *filename);
void save_png(const imgdta &image, std::ostream &os);
#endif

#ifdef WITH_JPEG
imgdta load_jpeg(const char *filename);
void save_jpeg(const imgdta &image, const char *filename, int quality = 80);
#endif

#ifdef WITH_TIFF
imgdta load_tiff(const char *filename);
void save_tiff(const imgdta &image, const char *filename);
#endif

#ifdef WITH_EXR
imgdta load_exr(const char *filename);
void save_exr(const imgdta &image, const char *filename);
#endif

imgdta load_pbm(const char *filename);
void save_pbm(const imgdta &image, const char *filename);

imgdta load(const char *filename);
void save(const imgdta &image, const char *filename);

}

namespace image_manip {

image_f grayscale(const image_b &src);
template <typename T = uint8_t> image<T> rescale_half_gaussian(const image<T> &img, float sigma = 0.866025403784439f);
template <typename T = uint8_t> image<T> blur_gaussian(const image<T> &img, float sigma);
template <typename T = uint8_t> image<T> subtract(const image<T> &a, const image<T> &b);

}



// implementations
namespace image_manip {

template <typename T>
image<T> rescale_half_gaussian(const image<T> &img, float sigma) // from MVE
{
	int iw = img.width();
	int ih = img.height();
	int ic = img.channels();
	int ow = (iw + 1) >> 1;
	int oh = (ih + 1) >> 1;

	if (iw < 2 || ih < 2)
		throw std::invalid_argument("Invalid input image");

	image<T> out(ow, oh, ic);

	/*
	 * Weights w1 (4 center px), w2 (8 skewed px) and w3 (4 corner px).
	 * Weights can be normalized by dividing with (4*w1 + 8*w2 + 4*w3).
	 * Since the accumulator is used, normalization is not explicit.
	 */
	float w1 = std::exp(-0.5f / (2.0f * sigma * sigma));
	float w2 = std::exp(-2.5f / (2.0f * sigma * sigma));
	float w3 = std::exp(-4.5f / (2.0f * sigma * sigma));
	float wsum = (w1 + w2 + w3) * 4;

	int outpos = 0;
	int rowstride = iw * ic;
	for (int y = 0; y < oh; ++y) {
		/* Init the four row pointers. */
		int y2 = (int)y << 1;
		T const* row[4];
		row[0] = img.data(std::max(0, y2 - 1) * rowstride);
		row[1] = img.data(y2 * rowstride);
		row[2] = img.data(std::min((int)ih - 1, y2 + 1) * rowstride);
		row[3] = img.data(std::min((int)ih - 1, y2 + 2) * rowstride);

		for (int x = 0; x < ow; ++x) {
			/* Init four pixel positions for each row. */
			int x2 = (int)x << 1;
			int xi[4];
			xi[0] = std::max(0, x2 - 1) * ic;
			xi[1] = x2 * ic;
			xi[2] = std::min((int)iw - 1, x2 + 1) * ic;
			xi[3] = std::min((int)iw - 1, x2 + 2) * ic;

			for (int c = 0; c < ic; ++c) {
				float accum = 0;
				accum += row[0][xi[0] + c] * w3;
				accum += row[0][xi[1] + c] * w2;
				accum += row[0][xi[2] + c] * w2;
				accum += row[0][xi[3] + c] * w3;

				accum += row[1][xi[0] + c] * w2;
				accum += row[1][xi[1] + c] * w1;
				accum += row[1][xi[2] + c] * w1;
				accum += row[1][xi[3] + c] * w2;

				accum += row[2][xi[0] + c] * w2;
				accum += row[2][xi[1] + c] * w1;
				accum += row[2][xi[2] + c] * w1;
				accum += row[2][xi[3] + c] * w2;

				accum += row[3][xi[0] + c] * w3;
				accum += row[3][xi[1] + c] * w2;
				accum += row[3][xi[2] + c] * w2;
				accum += row[3][xi[3] + c] * w3;

				out[outpos++] = std::min((float)std::numeric_limits<T>::max(), accum / wsum);
			}
		}
	}

	return out;
}

template <typename T>
image<T> blur_gaussian(const image<T> &img, float sigma)
{
	int ks = std::ceil(sigma * 2.884f); // Cap kernel at 1/128
	std::vector<float> kernel(ks + 1);
	float ksum = 0;
	for (int i = 0; i <= ks; ++i) {
		ksum += kernel[i] = std::exp(-(((float)i * (float)i) / (T(2) * sigma * sigma)));
		if (i) ksum += kernel[i];
	}
	for (int i = 0; i <= ks; ++i) {
		kernel[i] /= ksum;
	}

	image<T> tmp(img.height(), img.width(), img.channels());
	for (int y = 0; y < img.height(); ++y) {
		for (int x = 0; x < img.width(); ++x) {
			for (int c = 0; c < img.channels(); ++c) {
				float val = 0;
				for (int i = -ks; i <= ks; ++i) {
					int xx = std::min(std::max(x + i, 0), (int)img.width() - 1), yy = y;
					val += img(xx, yy, c) * kernel[std::abs(i)];
				}
				tmp(y, x, c) = val; // write transposed
			}
		}
	}

	image<T> out(img.width(), img.height(), img.channels());
	for (int x = 0; x < img.width(); ++x) {
		for (int y = 0; y < img.height(); ++y) {
			for (int c = 0; c < img.channels(); ++c) {
				float val = 0;
				for (int i = -ks; i <= ks; ++i) {
					int yy = std::min(std::max(y + i, 0), (int)img.height() - 1), xx = x;
					val += tmp(yy, xx, c) * kernel[std::abs(i)]; // read transposed (better cache utilization)
				}
				out(x, y, c) = val;
			}
		}
	}

	return out;
}

template <typename T>
image<T> subtract(const image<T> &a, const image<T> &b)
{
	image<T> tmp(a.width(), a.height(), a.channels());
	image<T> out(a.width(), a.height(), a.channels());
	for (int y = 0; y < a.height(); ++y) {
		for (int x = 0; x < a.width(); ++x) {
			for (int c = 0; c < a.channels(); ++c) {
				out(x, y, c) = a(x, y, c) - b(x, y, c);
			}
		}
	}
	return out;
}

}