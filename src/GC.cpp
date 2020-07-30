#include "definitions.h"
#include "utils.h"
#include "PRNG.h"
#include "blake3.h"
#include "OT.h"
#include "GC.h"

namespace GC
{
	uint32_t Gate::curId;
	Label Gate::delta;

	Label Label::operator^(const Label& label)
	{
		Label ret;
		for (int i = 0; i < bufflen; i++)
			ret.buffer[i] = this->buffer[i] ^ label.buffer[i];
		return ret;
	}

	void Label::operator^=(const Label& label)
	{
		for (int i = 0; i < bufflen; i++)
			this->buffer[i] ^= label.buffer[i];
	}

	bool Label::operator==(const Label& label)
	{
		for (int i = 0; i < bufflen; i++)
			if (buffer[i] != label.buffer[i])
				return false;
		return true;
	}

	Label Label::random(bool isDelta)
	{
		Label ret;
		for (int i = 0; i < bufflen; i++)
			ret.buffer[i] = gRNG.nextUInt64();
		if (isDelta)
			ret.buffer[0] |= 1;
		return ret;
	}
	Label Label::H(const Label& a, uint32_t index)
	{
		Label ret;
		blake3_hasher hasher;
		blake3_hasher_init(&hasher);
		blake3_hasher_update(&hasher, a.buffer, kappa / 8);
		blake3_hasher_update(&hasher, &index, 4);
		blake3_hasher_finalize(&hasher, (uint8_t*)ret.buffer, kappa / 8);
		return ret;
	}
	Label Label::H(const Label& a, const Label& b, uint32_t index)
	{
		Label ret;
		blake3_hasher hasher;
		blake3_hasher_init(&hasher);
		blake3_hasher_update(&hasher, a.buffer, kappa / 8);
		blake3_hasher_update(&hasher, b.buffer, kappa / 8);
		blake3_hasher_update(&hasher, &index, 4);
		blake3_hasher_finalize(&hasher, (uint8_t*)ret.buffer, kappa / 8);
		return ret;
	}
	void InGate::generate()
	{
		if (gOwner)
		{
			genOutLabel = Label::random();
			addComm(kappa); // sends the label
			if (input)
				evaOutLabel = genOutLabel ^ delta;
			else
				evaOutLabel = genOutLabel;
		}
		else
			OT::correlatedOT(input, delta.getBuffPt(), kappa / 8, genOutLabel.getBuffPt(), evaOutLabel.getBuffPt());
	}
	void GHAGate::generate(bool va)
	{
		Label hb0 = Label::H(b->genOutLabel, id);
		Label hb1 = Label::H(b->genOutLabel ^ delta, id);
		if (va)
			hb1 ^= delta;
		genOutLabel = b->genOutLabel.getPoint() ? hb1 : hb0;
		msg = hb1 ^ hb0;
	}
	void GHAGate::evaluate()
	{
		// The evaluator receives msg
		addComm(kappa);
		evaOutLabel = Label::H(b->evaOutLabel, id);
		if (b->evaOutLabel.getPoint()) // permute
			evaOutLabel ^= msg ;
	}

	void EHAGate::generate()
	{
		Label a0 = a->genOutLabel;
		Label a1 = a0 ^ delta;
		if (forHalfAND && a->genOutLabel.getPoint())
			std::swap(a0, a1);
		genOutLabel = Label::H(a0, id);
		msg = Label::H(a1, id) ^ genOutLabel ^ b->genOutLabel;
	}
	void EHAGate::evaluate()
	{
		// The evaluator receives msg
		addComm(kappa);
		bool input = forHalfAND ? a->evaOutLabel.getPoint() : dynamic_cast<InGate*> (a)->input;
		evaOutLabel = Label::H(a->evaOutLabel, id);
		if (input)
			evaOutLabel ^= msg ^ b->evaOutLabel;
	}

	ANDGate::ANDGate(Gate* a, Gate* b): Gate(a, b)
	{
		gha = new GHAGate(a);
		eha = new EHAGate(b, a, true);
		xorg = new XORGate(gha, eha);
	}
	void ANDGate::generate()
	{
		gha->generate(b->genOutLabel.getPoint());
		eha->generate();
		xorg->generate();
		genOutLabel = xorg->genOutLabel;

	}
	void ANDGate::evaluate()
	{
		gha->evaluate();
		eha->evaluate();
		xorg->evaluate();
		evaOutLabel = xorg->evaOutLabel;
		if (evaOutLabel != genOutLabel && evaOutLabel != (genOutLabel ^ delta))
		{
			throw "Error!";
		}
	}
	ANDGate::~ANDGate()
	{
		delete gha;
		delete eha;
		delete xorg;
	}

	void GC(const std::vector<Gate*>& gates, const std::vector<Gate*>& outGates, std::vector<bool>& z1, std::vector<bool>& z2)
	{
		Gate::curId = 0;
		Gate::delta = Label::random(true); // generate delta
		for (auto gate : gates) // One round communication 
			gate->generate();
		for (auto gate : gates) // One round communication
			gate->evaluate();
		size_t outSize = outGates.size();
		z1.resize(outSize);
		z2.resize(outSize);
		for (size_t i = 0; i < outSize; i++)
		{
			z1[i] = outGates[i]->genOutLabel.getPoint();
			z2[i] = outGates[i]->evaOutLabel.getPoint();
		}
	}
}