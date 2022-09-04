#include "image.h"

#define IMG_LOOP(E) for (std::size_t i = 0; i < src.pixels(); ++i) { dst[i] = E; }

namespace image_manip {

image_f grayscale(const image_b &src)
{
	image_f dst(src.width(), src.height());

	switch (src.channels()) {
	case 1:
		IMG_LOOP(src[i] / 255.f);
		break;
	case 2:
		IMG_LOOP((src.at(i, 0) / 255.f) * (src.at(i, 1) / 255.f));
		break;
	case 3:
		IMG_LOOP(src.at(i, 0) / 255.f * .3f + src.at(i, 1) / 255.f * .59f + src.at(i, 2) / 255.f * .11f);
		break;
	case 4:
		IMG_LOOP((src.at(i, 0) / 255.f * .3f + src.at(i, 1) / 255.f * .59f + src.at(i, 2) / 255.f * .11f) * (src.at(i, 1) / 255.f));
		break;
	}

	return dst;
}

}