#include <vector>

namespace GC
{
	class Label // A wire label of a gate
	{
	public:
		Label operator^(const Label& label);
		void operator^=(const Label& label);
		bool operator==(const Label& label);
		bool operator!=(const Label& label) { return !((*this) == label); }
		static Label random(bool isDelta = false);
		static Label H(const Label& a, uint32_t index);
		static Label H(const Label& a, const Label& b, uint32_t index);
		inline bool getPoint() { return buffer[0] & 1; } // interprete a bit as point and permute bit
		inline uint8_t* getBuffPt() { return (uint8_t*)buffer; }
	private:
		static const int bufflen = kappa / 64;
		uint64_t buffer[bufflen];
	};

	class Gate
	{
	public:
		static uint32_t curId; // current id to be allocated
		static Label delta;
		Gate(Gate* a, Gate* b) : a(a), b(b), id(curId++) {}
		Label genOutLabel; // Generator owns it
		Label evaOutLabel; // Evaluator owns it
		virtual void generate() = 0; // generator generates labels
		virtual void evaluate() = 0; // evaluator evaluates labels
		//virtual ~Gate() = 0;
	protected:
		uint32_t id;
		Gate* a, * b;
	};

	class InGate : public Gate // The input gate with input only a bit
	{
	private:
		bool gOwner; // Is generator the owner of the input?
	public:
		InGate(bool input, bool gOwner) : Gate(NULL,NULL), input(input), gOwner(gOwner) {}
		bool input; // the input bit
		void generate();
		void evaluate() {}; // This work has been done during generate phase
	};

	class XORGate : public Gate // a free XOR gate
	{
	public:
		XORGate(Gate* a, Gate* b) :Gate(a, b) {}
		void generate() { genOutLabel = a->genOutLabel ^ b->genOutLabel; }
		void evaluate() { evaOutLabel = a->evaOutLabel ^ b->evaOutLabel; }
	};

	class GHXGate : public Gate // a free XOR gate with va know by garbler
	{
	private:
		bool a;
	public:
		GHXGate(bool a, Gate* b) : Gate(NULL, b), a(a) {}
		void generate() { genOutLabel = a ? (b->genOutLabel ^ delta) : b->genOutLabel; }
		void evaluate() { evaOutLabel = b->evaOutLabel; }
	};

	class GHAGate : public Gate // generator half gate
	{
	private:
		bool va;
		Label msg; // message to be sent from generator to evaluator 
	public:
		GHAGate(Gate* b) :Gate(NULL, b), va(false) {} // for Half AND
		GHAGate(bool va, Gate* b) :Gate(NULL, b), va(va) {}
		void generate() { generate(va); }
		void generate(bool va); // for Half AND
		void evaluate();
	};

	class EHAGate : public Gate // evaluator half gate
	{
	private:
		bool forHalfAND; // this setting is for half and only
		Label msg; // message to be sent from generator to evaluator 
	public:
		EHAGate(Gate* a, Gate* b) :Gate(a, b),forHalfAND(false) {}
		EHAGate(Gate* a, Gate* b, bool forHalfAND) : Gate(a, b), forHalfAND(forHalfAND) {}
		void generate();
		void evaluate();
	};

	class ANDGate : public Gate // A full AND gate
	{
	private:
		GHAGate* gha;
		EHAGate* eha;
		XORGate*xorg;
	public:
		ANDGate(Gate* a, Gate* b);
		void generate();
		void evaluate();
		~ANDGate();
	};

	// Note that labels are generated and evaluated in the given order
	// outGates is the subset of gates whose outLabels are outputs of GC
	// Outputs are Boolean shared
	void GC(const std::vector<Gate*>& gates, const std::vector<Gate*>& outGates, std::vector<bool>& z1, std::vector<bool>& z2);

}