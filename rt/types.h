#pragma once

#include <vector>
#include <unordered_map>

namespace types {

static const int SIZES[11] = { 4, 8, 8, 8, 4, 4, 2, 2, 1, 1, 0 };
enum Type { FLOAT, DOUBLE, ULONG, LONG, UINT, INT, USHORT, SHORT, UCHAR, CHAR, NONE };
enum Interp {
	PX, PY, PZ, PW,
	NX, NY, NZ, NW,
	CR, CG, CB, CA,
	CAR, CAG, CAB, CAA, CAC,
	CDR, CDG, CDB, CDA, CDC,
	CSR, CSG, CSB, CSA, CSP, CSC,
	TU, TV, TW,
	SCALE, CONFIDENCE,

	OTHER, };


}

struct VtxBuf {
	std::vector<types::Type> types;
	std::vector<uint32_t> offsets;
	std::vector<uint8_t> buf;
	uint32_t idx_interp[types::OTHER];

	int curoff;
	int i;

	VtxBuf() : offsets(1, 0), i(0)
	{
		for (int i = 0; i < types::OTHER; ++i) idx_interp[i] = -1u;
	}

	void describe(types::Type t, types::Interp i)
	{
		idx_interp[i] = types.size();
		describe(t);
	}
	void describe(types::Type t)
	{
		types.push_back(t);
		offsets.push_back(offsets.back() + types::SIZES[t]);
		buf.resize(offsets.back());
	}

	uint8_t *getdst()
	{
		return buf.data() + offsets[i++];
	}

	void reset()
	{
		i = 0;
	}

	bool has(types::Interp i) const
	{
		return idx_interp[i] != -1u;
	}

	template <typename T>
	static T extract(const uint8_t *buf)
	{
		return *(T*)buf;
	}

	template <typename T>
	T get(types::Interp i) const
	{
		int idx = idx_interp[i];
		const uint8_t *b = buf.data() + offsets[idx];
		switch (types[idx]) {
		case types::NONE:   return T();
		case types::CHAR:   return (T)extract<int8_t>  (b);
		case types::UCHAR:  return (T)extract<uint8_t> (b);
		case types::SHORT:  return (T)extract<int16_t> (b);
		case types::USHORT: return (T)extract<uint16_t>(b);
		case types::INT:    return (T)extract<int32_t> (b);
		case types::UINT:   return (T)extract<uint32_t>(b);
		case types::FLOAT:  return (T)extract<float>   (b);
		case types::DOUBLE: return (T)extract<double>  (b);
		}
	}
};


struct Face {
	uint32_t idx[3];

	Face()
	{}
	Face(uint32_t a, uint32_t b, uint32_t c) : idx{ a, b, c }
	{}
	Face(const Face &face) : Face(face.idx[0], face.idx[1], face.idx[2])
	{}
	Face(const volatile Face &face) : idx{ face.idx[0], face.idx[1], face.idx[2] }
	{}
// 	Face plusoff(uint32_t off) const
// 	{
// 		return Face{ idx[0] + off, idx[1] + off, idx[2] + off };
// 	}
};
struct Face16 {
	uint16_t idx[3];
};

struct Mesh {
	typedef std::vector<Face> Faces;
	Faces faces;
	enum VtxDataType { VTX, NRM, TEX };
	std::unordered_map<uint32_t, std::vector<uint32_t>> fan_faces;
	void add_face(uint32_t a, uint32_t b, uint32_t c)
	{
		uint32_t x[] = { a, b, c };
		for (int z = 0; z < 3; ++z) {
			const std::vector<uint32_t> &cur_faces = fan_faces[x[z]];
			for (int i = 0; i < cur_faces.size(); ++i) {
				for (int j = 0; j < 3; ++j) {
					if (a == faces[cur_faces[i]].idx[j % 3] && b == faces[cur_faces[i]].idx[(j + 1) % 3] && c == faces[cur_faces[i]].idx[(j + 2) % 3]) {
						return;
					}
				}
			}
		}
		fan_faces[a].push_back(faces.size());
		fan_faces[b].push_back(faces.size());
		fan_faces[c].push_back(faces.size());
		faces.push_back(Face{ a, b, c });
	}
	void clear_tmp()
	{
		fan_faces = std::unordered_map<uint32_t, std::vector<uint32_t>>();
	}
	virtual void add_vtx(const VtxBuf &buf, VtxDataType t = VTX) = 0;
	virtual void validate(const VtxBuf &buf) const = 0;
	virtual std::size_t size_vertices() const = 0;
};
