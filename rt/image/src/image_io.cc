#include "image.h"

#include <stdexcept>
#include <fstream>
#include <string>
#include <bitset>
#include <limits>
#include <cstdio>

#ifdef WITH_PNG
#include <png.h>
#endif
#ifdef WITH_JPEG
#include <jpeglib.h>
#endif
#ifdef WITH_TIFF
#include <tiff.h>
#include <tiffio.h>
#endif
#ifdef WITH_EXR
#include <ImfOutputFile.h>
#include <ImfInputFile.h>
#include <ImfChannelList.h>
#include <ImfArray.h>
#endif

#undef max

namespace image_io {

#ifdef WITH_PNG
imgdta load_png(const char *filename)
{
	FILE *fp = std::fopen(filename, "rb");
	if (!fp) throw std::runtime_error("Cannot open PNG file");

	png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
	if (!png) {
		std::fclose(fp);
		throw std::runtime_error("Could not read PNG version");
	}

	png_infop info = png_create_info_struct(png);
	if (!info) {
		png_destroy_read_struct(&png, &info, nullptr);
		std::fclose(fp);
		throw std::runtime_error("Could not read PNG info");
	}

	if (setjmp(png_jmpbuf(png))) {
		png_destroy_read_struct(&png, &info, nullptr);
		std::fclose(fp);
		throw std::runtime_error("Internal libpng error");
	}

	png_init_io(png, fp);
	png_read_info(png, info);

	png_byte color_type = png_get_color_type(png, info);
	png_byte bit_depth = png_get_bit_depth(png, info);

	imgdta image(png_get_image_width(png, info), png_get_image_height(png, info), png_get_channels(png, info), bit_depth / 8);

// 	if (bit_depth == 16) png_set_strip_16(png);
	if (color_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png);
	if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) png_set_expand_gray_1_2_4_to_8(png);
	if (png_get_valid(png, info, PNG_INFO_tRNS)) png_set_tRNS_to_alpha(png);

	png_read_update_info(png, info);

	std::vector<png_bytep> row_pointers(image.height());
	for (std::size_t y = 0; y < image.height(); ++y) {
		row_pointers[y] = (png_byte*)image.ptr2d(0, y);
	}
	png_read_image(png, row_pointers.data());
	if (image.bps() == 2) {
		for (int i = 0; i < image.size(); ++i) {
			std::swap(image.data()[i * 2], image.data()[i * 2 + 1]);
		}
	} else if (image.bps() == 4) {
		for (int i = 0; i < image.size(); ++i) {
			std::swap(image.data()[i * 2], image.data()[i * 2 + 3]);
			std::swap(image.data()[i * 2 + 1], image.data()[i * 2 + 2]);
		}
	}

	png_destroy_read_struct(&png, &info, nullptr);
	std::fclose(fp);

	return image;
}

void save_png(const imgdta &image, const char *filename) // TODO: PNG16
{
	FILE *fp = std::fopen(filename, "wb");
	if (!fp) throw std::runtime_error("Cannot open PNG file");

	png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
	if (!png) {
		std::fclose(fp);
		throw std::runtime_error("Cannot create PNG write struct");
	}

	png_infop info = png_create_info_struct(png);
	if (!info) {
		png_destroy_write_struct(&png, &info);
		std::fclose(fp);
		throw std::runtime_error("Cannot create PNG info");
	}

	if (setjmp(png_jmpbuf(png))) {
		png_destroy_write_struct(&png, &info);
		std::fclose(fp);
		throw std::runtime_error("Internal libpng error");
	}

	png_init_io(png, fp);

	png_byte color_type;
	switch (image.channels()) {
		case 1: color_type = PNG_COLOR_TYPE_GRAY; break;
		case 2: color_type = PNG_COLOR_TYPE_GRAY_ALPHA; break;
		case 3: color_type = PNG_COLOR_TYPE_RGB; break;
		case 4: color_type = PNG_COLOR_TYPE_RGB_ALPHA; break;
		default: throw std::runtime_error("Cannot determine image color type");
	}
	png_set_IHDR(png, info, image.width(), image.height(), image.bps() * 8, color_type, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	png_write_info(png, info);

	std::vector<png_bytep> row_pointers(image.height());
	std::vector<uint8_t> data(image.bytes());
	if (image.bps() == 2)      for (int i = 0; i < image.size(); ++i) png_save_uint_16((png_bytep)((uint16_t*)data.data() + i), ((uint16_t*)image.data())[i]);
	else if (image.bps() == 4) for (int i = 0; i < image.size(); ++i) png_save_uint_32((png_bytep)((uint32_t*)data.data() + i), ((uint32_t*)image.data())[i]);
	else                       for (int i = 0; i < image.size(); ++i) data[i] = image.data()[i];
	for (std::size_t y = 0; y < image.height(); ++y) {
		row_pointers[y] = (png_byte*)(data.data() + y * image.pitch());
	}

	png_write_image(png, row_pointers.data());
	png_write_end(png, nullptr);

	png_destroy_write_struct(&png, &info);
	std::fclose(fp);
}
void png_write_data(png_structp pngPtr, png_bytep data, png_size_t length) {
	png_voidp a = png_get_io_ptr(pngPtr);
	((std::ostream*)a)->write((const char*)data, length);
}
void png_flush(png_structp pngPtr) {
	png_voidp a = png_get_io_ptr(pngPtr);
	((std::ostream*)a)->flush();
}
void save_png(const imgdta &image, std::ostream &os) // TODO: PNG16
{
	png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
	if (!png) {
// 		std::fclose(fp);
		throw std::runtime_error("Cannot create PNG write struct");
	}

	png_infop info = png_create_info_struct(png);
	if (!info) {
		png_destroy_write_struct(&png, &info);
// 		std::fclose(fp);
		throw std::runtime_error("Cannot create PNG info");
	}

	if (setjmp(png_jmpbuf(png))) {
		png_destroy_write_struct(&png, &info);
// 		std::fclose(fp);
		throw std::runtime_error("Internal libpng error");
	}

// 	png_init_io(png, fp);
	png_set_write_fn(png, (png_voidp)&os, png_write_data, png_flush);

	png_byte color_type;
	switch (image.channels()) {
		case 1: color_type = PNG_COLOR_TYPE_GRAY; break;
		case 2: color_type = PNG_COLOR_TYPE_GRAY_ALPHA; break;
		case 3: color_type = PNG_COLOR_TYPE_RGB; break;
		case 4: color_type = PNG_COLOR_TYPE_RGB_ALPHA; break;
		default: throw std::runtime_error("Cannot determine image color type");
	}
	png_set_IHDR(png, info, image.width(), image.height(), image.bps() * 8, color_type, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	png_write_info(png, info);

	std::vector<png_bytep> row_pointers(image.height());
	std::vector<uint8_t> data(image.bytes());
	if (image.bps() == 2)      for (int i = 0; i < image.size(); ++i) png_save_uint_16((png_bytep)((uint16_t*)data.data() + i), ((uint16_t*)image.data())[i]);
	else if (image.bps() == 4) for (int i = 0; i < image.size(); ++i) png_save_uint_32((png_bytep)((uint32_t*)data.data() + i), ((uint32_t*)image.data())[i]);
	else                       for (int i = 0; i < image.size(); ++i) data[i] = image.data()[i];
	for (std::size_t y = 0; y < image.height(); ++y) {
		row_pointers[y] = (png_byte*)(data.data() + y * image.pitch());
	}

	png_write_image(png, row_pointers.data());
	png_write_end(png, nullptr);

	png_destroy_write_struct(&png, &info);
// 	std::fclose(fp);
}
#endif

#ifdef WITH_JPEG
void jpg_error_handler(j_common_ptr)
{
	throw std::runtime_error("JPEG format not recognized");
}
void jpg_message_handler(j_common_ptr, int msg_level)
{
	if (msg_level < 0) throw std::runtime_error("JPEG data corrupt");
}

imgdta load_jpeg(const char *filename)
{
	FILE *fp = std::fopen(filename, "rb");
	if (!fp) throw std::runtime_error("Cannot open JPEG file");

	jpeg_decompress_struct cinfo;
	jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jerr.error_exit = &jpg_error_handler;
	jerr.emit_message = &jpg_message_handler;

	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, fp);
	jpeg_read_header(&cinfo, TRUE);

	imgdta image(cinfo.image_width, cinfo.image_height, cinfo.out_color_space == JCS_RGB ? 3 : 1);

	jpeg_start_decompress(&cinfo);

	unsigned char *data = image.data();
	while (cinfo.output_scanline < cinfo.output_height) {
		jpeg_read_scanlines(&cinfo, &data, 1);
		data += image.channels() * cinfo.output_width;
	}

	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	std::fclose(fp);

	return image;
}

void save_jpeg(const imgdta &image, const char *filename, int quality)
{
	if (image.channels() != 1 && image.channels() != 3) throw std::runtime_error("Invalid number of channels");

	FILE *fp = std::fopen(filename, "wb");
	if (!fp) throw std::runtime_error("Cannot open JPEG file");

	jpeg_compress_struct cinfo;
	jpeg_error_mgr jerr;

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, fp);

	cinfo.image_width = image.width();
	cinfo.image_height = image.height();
	cinfo.input_components = image.channels();
	cinfo.in_color_space = image.channels() == 1 ? JCS_GRAYSCALE : JCS_RGB;

	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality, TRUE);
	jpeg_start_compress(&cinfo, TRUE);

	while (cinfo.next_scanline < cinfo.image_height) {
		JSAMPROW row_pointer = (JSAMPROW)image.ptr2d(0, cinfo.next_scanline);
		jpeg_write_scanlines(&cinfo, &row_pointer, 1);
	}

	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
	std::fclose(fp);
}
#endif

#ifdef WITH_TIFF
void tiff_error_handler (char const*, char const *fmt, va_list ap)
{
	char msg[2048];
	::vsprintf(msg, fmt, ap);
	throw std::runtime_error(msg);
}


imgdta load_tiff(const char *filename)
{
	TIFFSetWarningHandler(nullptr);
	TIFFSetErrorHandler(tiff_error_handler);

	TIFF* tif = TIFFOpen(filename, "r");
	if (!tif) throw std::runtime_error("TIFF file format not recognized");
	uint32_t width, height;
	uint16_t channels, bits;
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
	TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &channels);
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits);

	imgdta image(width, height, channels, bits / 8);
	uint32_t rowstride = TIFFScanlineSize(tif);
	for (uint32_t row = 0; row < height; row++) {
		tdata_t row_pointer = image.data() + row * rowstride;
		TIFFReadScanline(tif, row_pointer, row);
	}
	TIFFClose(tif);

	return image;
}

void save_tiff(const imgdta &image, const char *filename)
{
	TIFF* tif = TIFFOpen(filename, "w");
	if (!tif) throw std::runtime_error("Unknown TIFF file error");

	uint32_t width = image.width();
	uint32_t height = image.height();
	uint32_t channels = image.channels();
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, channels);
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, image.bps() * 8);
// 	TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
	TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

	int64_t ret = TIFFWriteEncodedStrip(tif, 0, (uint8_t*)image.data(), image.bytes());

	TIFFClose(tif);

	if (ret < 0) throw std::runtime_error("Error writing TIFF image");
}
#endif

#ifdef WITH_EXR
const char *exr_channel_names[4][4] = { { "Y", "", "", "" }, { "Y", "A", "", "" }, { "R", "G", "B", "" }, { "R", "G", "B", "A" } };

imgdta load_exr(const char *filename)
{
	Imf::InputFile file(filename);
	const Imf::Header &header = file.header();
	enum Channels { R = 1, G = 2, B = 4, A = 8, Y = 16 };
	int channels = 0;
	for (Imf::ChannelList::ConstIterator it = header.channels().begin(); it != header.channels().end(); ++it) {
		switch (it.name()[0]) {
		case 'R': channels |= R; break;
		case 'G': channels |= G; break;
		case 'B': channels |= B; break;
		case 'A': channels |= A; break;
		case 'Y': channels |= Y; break;
		default: throw std::runtime_error("Unsupported channel");
		}
	}
	std::size_t nc = 0;
	if (channels & R && channels & G && channels & B && channels & A) nc = 4;
	else if (channels & R && channels & G && channels & B) nc = 3;
	else if (channels & Y && channels & A) nc = 2;
	else if (channels & Y) nc = 1;
	else throw std::runtime_error("Unsupported channels");

	Imath::Box2i dw = header.dataWindow();
	imgdta image(dw.max.x - dw.min.x + 1, dw.max.y - dw.min.y + 1, nc, 4, true);

	Imf::FrameBuffer fb;
	for (std::size_t c = 0; c < nc; ++c) {
		fb.insert(exr_channel_names[nc][c], Imf::Slice(Imf::FLOAT, (char*)((char*)image.data() + sizeof(float) * c), sizeof(float) * nc, sizeof(float) * image.width() * nc, 1, 1, 0.0));
	}

	file.setFrameBuffer(fb);
	file.readPixels(dw.min.y, dw.max.y);

	return image;
}

void save_exr(const imgdta &image, const char *filename)
{
	if (image.channels() < 1 || image.channels() > 4 || !image.is_float()) throw std::runtime_error("Invalid number of channels");

	Imf::FrameBuffer fb;
	Imf::Header header(image.width(), image.height());
	std::size_t nc = image.channels();
	for (std::size_t c = 0; c < nc; ++c) {
		fb.insert(exr_channel_names[nc][c], Imf::Slice(Imf::FLOAT, (char*)((char*)image.data() + sizeof(float) * c), sizeof(float) * nc, sizeof(float) * image.width() * nc, 1, 1, 0.0));
		header.channels().insert(exr_channel_names[nc][c], Imf::Channel(Imf::FLOAT));
	}

	Imf::OutputFile outFile(filename, header);
	outFile.setFrameBuffer(fb);
	outFile.writePixels(image.height());
}
#endif

imgdta load_pbm(const char *filename)
{
	std::ifstream is(filename, std::ios_base::binary);
	std::string magic;
	is >> magic;
	std::size_t channels;
	enum Header { FLOAT = 1, BIT = 2, ASCII = 4, RGB = 8 };
	int header = 0;
	if (magic == "Pf") {
		header |= FLOAT;
	} else if (magic == "PF") {
		header |= FLOAT;
		header |= RGB;
	} else if (magic == "P1") {
		header |= BIT;
		header |= ASCII;
	} else if (magic == "P2") {
		header |= ASCII;
	} else if (magic == "P3") {
		header |= RGB;
		header |= ASCII;
	} else if (magic == "P4") {
		header |= BIT;
	} else if (magic == "P5") {
	} else if (magic == "P6") {
		header |= RGB;
	}

	float scale = 1, f_val;
	std::size_t i_val, w, h;
	is >> w >> h;
	if (!(header & BIT)) is >> scale;
	scale = 1. / scale;
	if (!(header & FLOAT)) scale *= 255;

	imgdta image(w, h, header & RGB ? 3 : 1, (header & FLOAT) ? 4 : 1, header & FLOAT);

	is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	if (!(header & ASCII)) {
		if (header & BIT) {
			uint8_t *cur = image.data();
			std::vector<std::bitset<8>> row((image.width() + 7) / 8);
			for (std::size_t y = 0; y < image.height(); ++y) {
				is.read((char*)row.data(), row.size());
				for (std::size_t x = 0; x < image.width(); ++x) {
					*cur++ = !row[x >> 3][7 - (x & 7)];
				}
			}
		} else {
			is.read((char*)image.data(), image.bytes());
		}
	} else if (header & FLOAT) {
		float *cur = (float*)image.data();
		for (std::size_t i = 0; i < image.size(); ++i) {
			is >> *cur++;
		}
	} else {
		uint8_t *cur = image.data();
		for (std::size_t i = 0; i < image.size(); ++i) {
			is >> i_val;
			*cur++ = header & BIT ? !i_val : i_val;
		}
	}

	if (!(header & FLOAT)) {
		uint8_t *cur = image.data();
		for (std::size_t i = 0; i < image.size(); ++i) {
			*cur++ *= scale;
		}
	}

	return image;
}

void save_pbm(const imgdta &image, const char *filename)
{
	std::ofstream os(filename, std::ios_base::binary);
	if (image.is_float() && image.bps() == 4 && image.channels() == 1) os << "Pf\n";
	else if (image.is_float() && image.bps() == 4 && image.channels() == 3) os << "PF\n";
	else if (image.channels() == 1) os << "P5\n";
	else if (image.channels() == 3) os << "P6\n";
	else throw std::runtime_error("Unsupported number of channels");

	os << image.width() << " " << image.height() << "\n";
	if (image.is_float()) os << "1\n";
	else os << "255\n";
	os.write((const char*)image.data(), image.bytes());
}

#define MVEI_FILE_SIGNATURE "\211MVE_IMAGE\n"
#define MVEI_FILE_SIGNATURE_LEN 11
#define MVEI_MAX_PIXEL_AMOUNT (16384 * 16384) /* 2^28 */

imgdta load_mvei(const char *filename)
{
	std::ifstream is(filename, std::ios_base::binary);
	char sig[MVEI_FILE_SIGNATURE_LEN];
	is.read(sig, MVEI_FILE_SIGNATURE_LEN);
	if (!std::equal(sig, sig + MVEI_FILE_SIGNATURE_LEN, MVEI_FILE_SIGNATURE)) throw std::runtime_error("Invalid MVEI signature");
	int32_t w, h, nc, t;
	is.read((char*)&w, sizeof(int32_t));
	is.read((char*)&h, sizeof(int32_t));
	is.read((char*)&nc, sizeof(int32_t));
	is.read((char*)&t, sizeof(int32_t));
	if (t < 1 || t > 10) throw std::runtime_error("Invalid MVEI type");

	static const int bytes[] = { 0, 1, 2, 4, 8, 1, 2, 4, 8, 4, 8 };
	imgdta image(w, h, nc, bytes[t], t == 9 || t == 10);

	is.read((char*)image.data(), image.size());

	return image;
}

void save_mvei(const imgdta &image, const char *filename)
{
	std::ofstream os(filename, std::ios_base::binary);
	os.write(MVEI_FILE_SIGNATURE, MVEI_FILE_SIGNATURE_LEN);

	static const int bps_i[] = { 0, 1, 2, 0, 3, 0, 0, 0, 4 }; // dont write signed
	static const int bps_f[] = { 0, 0, 0, 0, 9, 0, 0, 0, 10 };
	int32_t w = image.width(), h = image.height(), nc = image.channels(), t = image.is_float() ? bps_f[image.bps()] : bps_i[image.bps()];

	os.write((const char*)&w, sizeof(int32_t));
	os.write((const char*)&h, sizeof(int32_t));
	os.write((const char*)&nc, sizeof(int32_t));
	os.write((const char*)&t, sizeof(int32_t));

	os.write((const char*)image.data(), image.bytes());
}

enum Extension { PNG = 0x706e67, JPG = 0x6a7067, JPEG = 0x6a706567, TIF = 0x746966, TIFF = 0x74696666, EXR = 0x657872, PBM = 0x70626d, PGM = 0x70676d, PPM = 0x70706d, PFM = 0x70666d, MVEI = 0x6d766569 };

Extension get_extension(const char *filename)
{
	int dot = -1;
	for (int i = 0; filename[i]; ++i) {
		if (filename[i] == '.') dot = i;
	}
	uint32_t ext = 0;
	for (int i = dot + 1; filename[i]; ++i) {
		ext <<= 8;
		ext |= (uint8_t)::tolower(filename[i]);
	}
	return (Extension)ext;
}

imgdta load(const char *filename)
{
	switch (get_extension(filename)) {
#ifdef WITH_PNG
	case PNG:
		return load_png(filename);
#endif
#ifdef WITH_JPEG
	case JPG:
	case JPEG:
		return load_jpeg(filename);
#endif
#ifdef WITH_TIFF
	case TIF:
	case TIFF:
		return load_tiff(filename);
#endif
#ifdef WITH_EXR
	case EXR:
		return load_exr(filename);
#endif
	case PBM:
	case PGM:
	case PPM:
	case PFM:
		return load_pbm(filename);
	case MVEI:
		return load_mvei(filename);
	default:
		throw std::runtime_error("Unsupported file extension");
	}
}

void save(const imgdta &image, const char *filename)
{
	switch (get_extension(filename)) {
#ifdef WITH_PNG
	case PNG:
		return save_png(image, filename);
#endif
#ifdef WITH_JPEG
	case JPG:
	case JPEG:
		return save_jpeg(image, filename);
#endif
#ifdef WITH_TIFF
	case TIF:
	case TIFF:
		return save_tiff(image, filename);
#endif
#ifdef WITH_EXR
	case EXR:
		return save_exr(image, filename);
#endif
	case PBM:
	case PGM:
	case PPM:
	case PFM:
		return save_pbm(image, filename);
	case MVEI:
		return save_mvei(image, filename);
	default:
		throw std::runtime_error("Unsupported file extension");
	}
}

}