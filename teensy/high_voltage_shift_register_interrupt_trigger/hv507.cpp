#include "hv507.h"

void HV507::init(int8_t clockPin, int8_t latchPin, int8_t dataPin, int8_t polarityPin, int8_t blankingPin, int8_t directionPin)
{
  _clockPin = clockPin;
  _latchPin = latchPin;
  _dataPin = dataPin;
  _polarityPin = polarityPin;
  _blankingPin = blankingPin;
  _directionPin = directionPin;

  pinMode(_clockPin, OUTPUT);
  pinMode(_latchPin, OUTPUT);
  pinMode(_dataPin, OUTPUT);
  pinMode(_polarityPin, OUTPUT);
  pinMode(_blankingPin, OUTPUT);
  pinMode(_directionPin, OUTPUT);
  digitalWrite(_clockPin, LOW);  // don't assume what _clockPin is after each function, usu. HIGH though
  digitalWrite(_latchPin, LOW);  // _latchPin is LOW after every function
  digitalWrite(_dataPin, LOW);
  digitalWrite(_polarityPin, HIGH);  // _polarityPin is HIGH after each function
  digitalWrite(_blankingPin, HIGH);  // _blankingPin is HIGH after each function
  digitalWrite(_directionPin, LOW);  // LOW if dataPin connected to DIOA, HIGH if connected to DIOB
  
  _isInitialized = true;
  _invertedOutput = false;
}

void HV507::loadBitToRegister(uint8_t data)
{
	//digitalWrite(_latchPin, LOW);
  digitalWrite(_clockPin, LOW);
  digitalWrite(_dataPin, data);
  digitalWrite(_clockPin, HIGH); // Set clock pin low ready for a rising edge
}

void HV507::loadByteToRegister(uint8_t data)  // bits 7...0 = HVOUT 7...0
{
	if (_invertedOutput) { invertOutputs(); }
    for (int i = 0; i <= 7; i++)
    {
      digitalWrite(_clockPin, LOW);
		  digitalWrite(_dataPin, data & (1 << i));
      digitalWrite(_clockPin, HIGH);
    }
}

// data1 = HVOUT 0...15, data2 = HVOUT 16...31, data3 = HVOUT 32...45, data4 = HVOUT 46...63
void HV507::loadArrayToRegister(uint16_t data1, uint16_t data2, uint16_t data3, uint16_t data4)
{
	if (_invertedOutput) { invertOutputs(); }
  for (int i = 0; i <= 15; i++) {
    digitalWrite(_clockPin, LOW);
    if (data1 & (1 << i)) {
      digitalWrite(_dataPin, HIGH);
    } else {
      digitalWrite(_dataPin, LOW);
    }
    digitalWrite(_clockPin, HIGH);
  }
  for (int i = 0; i <= 15; i++) {
    digitalWrite(_clockPin, LOW);
    if (data2 & (1 << i)) {
      digitalWrite(_dataPin, HIGH);
    } else {
      digitalWrite(_dataPin, LOW);
    }
    digitalWrite(_clockPin, HIGH);
  }
  for (int i = 0; i <= 15; i++) {
    digitalWrite(_clockPin, LOW);
    if (data3 & (1 << i)) {
      digitalWrite(_dataPin, HIGH);
    } else {
      digitalWrite(_dataPin, LOW);
    }
    digitalWrite(_clockPin, HIGH);
  }
  for (int i = 0; i <= 15; i++) {
    digitalWrite(_clockPin, LOW);
    if (data4 & (1 << i)) {
      digitalWrite(_dataPin, HIGH);
    } else {
      digitalWrite(_dataPin, LOW);
    }
    digitalWrite(_clockPin, HIGH);
  }
}

// data1 = HVOUT 0...15, data2 = HVOUT 16...31, data3 = HVOUT 32...45, data4 = HVOUT 46...63
void HV507::loadArrayDirectlyToOutput(uint16_t data1, uint16_t data2, uint16_t data3, uint16_t data4)
{
  digitalWrite(_latchPin, HIGH);
	loadArrayToRegister(data1, data2, data3, data4);
	digitalWrite(_latchPin, LOW);
}

// load the values in the shift register to the high-voltage output
void HV507::loadInverseRegisterToOutput()
{
	// You'll see the output start to invert initially (because this matches invertOutputs())
	// On an Arduino Micro, the _latchPin commands starts ~5us after this, so you'll see 
	// the inverse of the output occur on the output ~10us after you call this function
	// (because the latches send their data on the LE pin's falling edge)
	digitalWrite(_polarityPin, LOW);
  digitalWrite(_latchPin, HIGH);
  digitalWrite(_latchPin, LOW);
	digitalWrite(_polarityPin, HIGH);
}

void HV507::loadRegisterToOutput()
{
	if (_invertedOutput) { invertOutputs(); }
    digitalWrite(_latchPin, HIGH);
    digitalWrite(_latchPin, LOW); // LE on for ~5.2us on an Arduino Micro (seems to work well)
}


void HV507::allOutputsHigh()
{
	if (_invertedOutput) { invertOutputs(); }
  digitalWrite(_polarityPin, LOW);
	digitalWrite(_blankingPin, LOW);
  digitalWrite(_blankingPin, HIGH);
	digitalWrite(_polarityPin, HIGH);
}

void HV507::allOutputsLow()
{
	if (_invertedOutput) { invertOutputs(); }
	digitalWrite(_blankingPin, LOW);
  digitalWrite(_blankingPin, HIGH);
}

void HV507::invertOutputs()
{
	if (_invertedOutput) {
		digitalWrite(_polarityPin, HIGH);
	} else {
		digitalWrite(_polarityPin, LOW);
	}
	_invertedOutput = !_invertedOutput;
}