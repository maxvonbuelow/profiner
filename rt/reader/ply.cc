/*
 * Copyright (C) 2017, Max von Buelow
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include "ply.h"

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <limits>
#include <numeric>
#include <algorithm>

#include "utils/io.h"
#include "utils/endian.h"

namespace reader {
namespace ply {

inline float uint2float(uint32_t i)
{
	return *((float*)&i);
}
inline double uint2double(uint64_t i)
{
	return *((double*)&i);
}
template <typename T>
T is_read_type(std::istream &is) {
	T val;
	is.read((char*)&val, sizeof(T));
	return val;
}

enum Fmt { ASCII, BIN_LE, BIN_BE };
typedef types::Type Type;

static std::unordered_map<std::string, int> lut_well_known = {
	{ "x", types::PX }, { "y", types::PY }, { "z", types::PZ }, { "w", types::PW },
	{ "nx", types::NX }, { "ny", types::NY }, { "nz", types::NZ }, { "nw", types::NW },
	{ "red", types::CR }, { "green", types::CG }, { "blue", types::CB }, { "alpha", types::CA },
	{ "ambient_red", types::CAR }, { "ambient_green", types::CAG }, { "ambient_blue", types::CAB }, { "ambient_alpha", types::CAA }, { "ambient_coeff", types::CAC },
	{ "diffuse_red", types::CDR }, { "diffuse_green", types::CDG }, { "diffuse_blue", types::CDB }, { "diffuse_alpha", types::CDA }, { "diffuse_coeff", types::CDC },
	{ "specular_red", types::CSR }, { "specular_green", types::CSG }, { "specular_blue", types::CSB }, { "specular_alpha", types::CSA }, { "specular_power", types::CSP }, { "specular_coeff", types::CSC },
	{ "u", types::TU }, { "tu", types::TU }, { "v", types::TV }, { "tv", types::TV }, { "tw", types::TW },
	{ "value", types::SCALE }, { "scale", types::SCALE }, { "confidence", types::CONFIDENCE },
};

inline Type str2type(const std::string &s)
{
	if (s == "float"  || s == "float32") return types::FLOAT;
	if (s == "double" || s == "float64") return types::DOUBLE;
	if (s == "uint"   || s == "uint32")  return types::UINT;
	if (s == "int"    || s == "int32")   return types::INT;
	if (s == "ushort" || s == "uint16")  return types::USHORT;
	if (s == "short"  || s == "int16")   return types::SHORT;
	if (s == "uchar"  || s == "uint8")   return types::UCHAR;
	if (s == "char"   || s == "int8")    return types::CHAR;
	throw std::runtime_error("Invalid data type");
};

struct Property {
	std::string name;
	Type type;
	Type list_len_type;

	inline Property(const std::string &_name, Type _type, Type _list_len_type = types::NONE) : name(_name), type(_type), list_len_type(_list_len_type)
	{}
};
struct Element : std::vector<Property> {
	std::string name;
	int len;
	std::unordered_map<std::string, int> propmap;

	inline Element(const std::string &_name, int _len) : name(_name), len(_len)
	{}

	inline int add(const std::string &name, Type type, Type list_len_type = types::NONE)
	{
		int idx = this->size();
		propmap[name] = idx;
		this->push_back(Property(name, type, list_len_type));
		return idx;
	}
	inline int operator[](const std::string &name) const
	{
		std::unordered_map<std::string, int>::const_iterator it = propmap.find(name);
		return it == propmap.end() ? -1 : it->second;
	}
	inline const Property &operator[](int i) const
	{
		return std::vector<Property>::operator[](i);
	}
	inline Property &operator[](int i)
	{
		return std::vector<Property>::operator[](i);
	}

	void init(VtxBuf &buf) const
	{
		for (int i = 0; i < this->size(); ++i) {
			const Property &prop = (*this)[i];
			std::unordered_map<std::string, int>::const_iterator it = lut_well_known.find(prop.name);
			if (prop.list_len_type != types::NONE) continue;
			if (it == lut_well_known.end())
				buf.describe(prop.type);
			else
				buf.describe(prop.type, (types::Interp)it->second);
		}
	}
};
struct Header : std::vector<Element> {
	Fmt fmt;
	std::unordered_map<std::string, int> elemmap;

	inline int add(const std::string &name, int len)
	{
		int idx = this->size();
		elemmap[name] = idx;
		this->push_back(Element(name, len));
		return idx;
	}

	inline int operator[](const std::string &name) const
	{
		std::unordered_map<std::string, int>::const_iterator it = elemmap.find(name);
		return it == elemmap.end() ? -1 : it->second;
	}
	inline const Element &operator[](int i) const
	{
		return std::vector<Element>::operator[](i);
	}
	inline Element &operator[](int i)
	{
		return std::vector<Element>::operator[](i);
	}
};

Header read_header(std::istream &is)
{
	std::string id;
	Header header;
	int curelem = -1;

	do {
		is >> id;
		if (id == "ply") {
			// #whatever
		} else if (id == "format") {
			std::string fmt;
			is >> fmt;
			if (fmt == "ascii") header.fmt = ASCII;
			else if (fmt == "binary_little_endian") header.fmt = BIN_LE;
			else if (fmt == "binary_big_endian") header.fmt = BIN_BE;
			else throw std::runtime_error("Invlaid format");
			util::skip_line(is);
		} else if (id == "comment") {
			util::skip_line(is);
		} else if (id == "element") {
			std::string name;
			is >> name;
			int len;
			is >> len;
			curelem = header.add(name, len);
		} else if (id == "property") {
			if (curelem == -1) throw std::runtime_error("Invlaid property");
			std::string type, name;
			types::Type list_len_type = types::NONE;

			is >> type;
			if (type == "list") {
				is >> type;
				list_len_type = str2type(type);
				is >> type;
			}
			is >> name;

			header[curelem].add(name, str2type(type), list_len_type);
		} else if (id != "end_header") {
			std::cout << "Skipping invalid PLY line " << id << std::endl;
			util::skip_line(is);
		}
	} while (id != "end_header");
	util::skip_line(is);

	return header;
}

template <typename T>
T cpval2ptr(unsigned char *ptr, T val)
{
	unsigned char *src = (unsigned char*)&val;
	std::copy(src, src + sizeof(T), ptr);
	return val;
}
struct ASCIIReader {
	uint64_t operator()(std::istream &is, unsigned char *dst, Type type)
	{
		int64_t val;
		uint64_t uval;
		double fval;
		switch (type) {
		case types::NONE:
			return 1;
		case types::CHAR:
			is >> val;
			return cpval2ptr<int8_t>(dst, val);
		case types::UCHAR:
			is >> uval;
			return cpval2ptr<uint8_t>(dst, uval);
		case types::SHORT:
			is >> val;
			return cpval2ptr<int16_t>(dst, val);
		case types::USHORT:
			is >> uval;
			return cpval2ptr<uint16_t>(dst, uval);
		case types::INT:
			is >> val;
			return cpval2ptr<int32_t>(dst, val);
		case types::UINT:
			is >> uval;
			return cpval2ptr<uint32_t>(dst, uval);
		case types::FLOAT:
			is >> fval;
			return cpval2ptr<float>(dst, fval);
		case types::DOUBLE:
			is >> fval;
			return cpval2ptr<double>(dst, fval);
		}
	}
};
struct BinBEReader {
	uint64_t operator()(std::istream &is, unsigned char *dst, Type type)
	{
		switch (type) {
		case types::NONE:   return 1;
		case types::CHAR:   return cpval2ptr<int8_t>  (dst,                     is_read_type<int8_t>(is));
		case types::UCHAR:  return cpval2ptr<uint8_t> (dst,                     is_read_type<uint8_t>(is));
		case types::SHORT:  return cpval2ptr<int16_t> (dst,             be16toh(is_read_type<int16_t>(is)));
		case types::USHORT: return cpval2ptr<uint16_t>(dst,             be16toh(is_read_type<uint16_t>(is)));
		case types::INT:    return cpval2ptr<int32_t> (dst,             be32toh(is_read_type<int32_t>(is)));
		case types::UINT:   return cpval2ptr<uint32_t>(dst,             be32toh(is_read_type<uint32_t>(is)));
		case types::FLOAT:  return cpval2ptr<float>   (dst, uint2float( be32toh(is_read_type<uint32_t>(is))));
		case types::DOUBLE: return cpval2ptr<double>  (dst, uint2double(be64toh(is_read_type<uint64_t>(is))));
		}
	}
};
struct BinLEReader {
	uint64_t operator()(std::istream &is, unsigned char *dst, Type type)
	{
		switch (type) {
		case types::NONE:   return 1;
		case types::CHAR:   return cpval2ptr<int8_t>  (dst,                     is_read_type<int8_t>(is));
		case types::UCHAR:  return cpval2ptr<uint8_t> (dst,                     is_read_type<uint8_t>(is));
		case types::SHORT:  return cpval2ptr<int16_t> (dst,             le16toh(is_read_type<int16_t>(is)));
		case types::USHORT: return cpval2ptr<uint16_t>(dst,             le16toh(is_read_type<uint16_t>(is)));
		case types::INT:    return cpval2ptr<int32_t> (dst,             le32toh(is_read_type<int32_t>(is)));
		case types::UINT:   return cpval2ptr<uint32_t>(dst,             le32toh(is_read_type<uint32_t>(is)));
		case types::FLOAT:  return cpval2ptr<float>   (dst, uint2float( le32toh(is_read_type<uint32_t>(is))));
		case types::DOUBLE: return cpval2ptr<double>  (dst, uint2double(le64toh(is_read_type<uint64_t>(is))));
		}
	}
};

template <typename R>
void readloop(std::istream &is, const Header &header, VtxBuf &buf, Mesh &mesh, R &&read)
{
	int face_idx = header["face"];
	int vtx_idx = header["vertex"];
	int vi_idx = header[face_idx]["vertex_indices"];
	uint32_t cnt = 0, cur = 0;
	for (int i = 0; i < header.size(); ++i) {
		if (i == face_idx || i == vtx_idx) cnt += header[i].len;
	}
	for (int i = 0; i < header.size(); ++i) {
		const Element &elem = header[i];
		uint8_t ign[8];

		if (i == face_idx || i == vtx_idx) {
			for (int j = 0; j < elem.len; ++j) {
				buf.reset();
				for (int k = 0; k < elem.size(); ++k) {
					const Property &prop = elem[k];
					if (prop.list_len_type != types::NONE) { // not an attribute
						uint64_t listlen = read(is, ign, prop.list_len_type);
						if (i == face_idx && k == vi_idx) { // ...but connectivity
							if (listlen > 0) {
								uint32_t first = read(is, ign, prop.type);
								if (listlen > 1) {
									uint32_t last = read(is, ign, prop.type);
									for (int l = 2; l < listlen; ++l) {
										uint32_t cur = read(is, ign, prop.type);
										mesh.add_face(first, last, cur);
										last = cur;
									}
								}
							}
						} else { // no connectivity and not an attribute: ignore
							for (int l = 0; l < listlen; ++l) {
								read(is, ign, prop.type);
							}
						}
					} else { // attribute
						if (i == face_idx) read(is, ign, prop.type);
						else read(is, buf.getdst(), prop.type);
					}
				}
				if (i == vtx_idx) mesh.add_vtx(buf);
			}
		} else {
			// ignore unknown elements
			for (int j = 0; j < elem.len; ++j) {
				for (int k = 0; k < elem.size(); ++k) {
					const Property &prop = elem[k];
					uint64_t listlen = read(is, ign, prop.list_len_type);
					for (int l = 0; l < listlen; ++l) {
						read(is, ign, prop.type);
					}
				}
			}
		}
	}
}

void read(std::istream &is, Mesh &mesh)
{
	Header header = read_header(is);

	VtxBuf buf;
	header[header["vertex"]].init(buf);
	mesh.validate(buf);

	// load mesh
	switch (header.fmt)
	{
	case ASCII:
		readloop(is, header, buf, mesh, ASCIIReader());
		break;
	case BIN_BE:
		readloop(is, header, buf, mesh, BinBEReader());
		break;
	case BIN_LE:
		readloop(is, header, buf, mesh, BinLEReader());
		break;
	}
}

}
}
