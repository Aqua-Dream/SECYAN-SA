#pragma once
#include <cstdint>

namespace OT
{
	void OT(uint8_t* msg0, uint8_t* msg1, int length, bool b, uint8_t* msgb);
	void GMW_OT(bool b, bool& msg0, bool& msg1, bool& msgb) ;
	void correlatedOT(bool b, const uint8_t* delta, int length, uint8_t* msg0, uint8_t* msgb) ;
	void shareMulCOT(bool b, uint16_t delta, uint16_t& msg0, uint16_t& msgb) ;
};