/***
Teensy 3.6 code for controlling the Microchip HV507 high voltage shift register

When the Teensy receives a trigger, it starts toggling the HV507 at a pre-defined bipolar
    square wave frequency. The code accepts two types of triggers: either the other Teensy
    triggers it by connecting the two Teensy's ground pins and linking the other Teensy's
    INTERRUPT_PIN to this Teensy's BUTTON_IN_PIN, or a push button can be connected between 
    this Teensy's BUTTON_IN_PIN and BUTTON_OUT_PIN (for other applications).
***/

#include "hv507.h"

const int PIN_CLK_EN = 24;
const int PIN_CLK = 25;
const int PIN_DIOB = 26;
const int PIN_DIR = 27;
const int PIN_BL = 28;
const int PIN_LE = 29;
const int PIN_POL = 30;
const int PIN_DIOA = 31;

const int LED_PIN = 13;
const int BUTTON_IN_PIN = 34;
const int BUTTON_OUT_PIN = 33;
const float FREQUENCY = 100;  // bipolar square wave drive frequency
const int DELAY = 1e6 / FREQUENCY / 2;
volatile int currIdx = 0;
bool forward = true;
volatile bool enabled = false;
volatile bool voltage_on_electrode1 = false;

HV507 hv507;
IntervalTimer timerBipolarSquareWave;

void setup() {
  pinMode(PIN_DIOB, INPUT);
  pinMode(BUTTON_OUT_PIN, OUTPUT);
  hv507.init(PIN_CLK, PIN_LE, PIN_DIOA, PIN_POL, PIN_BL, PIN_DIR);
  delay(100);
  hv507.loadArrayToRegister(0, 0, 0, 0);
  hv507.loadRegisterToOutput();
  pinMode(LED_PIN, OUTPUT);
  pinMode(BUTTON_IN_PIN, INPUT);
  digitalWrite(BUTTON_OUT_PIN, LOW);
  attachInterrupt(digitalPinToInterrupt(BUTTON_IN_PIN), hv507Toggle, CHANGE);
}

void loop() { }

void hv507Toggle() {
  if (digitalRead(BUTTON_IN_PIN) == HIGH) {
    timerBipolarSquareWave.end();
    digitalWrite(LED_PIN, HIGH);
    hv507.loadArrayToRegister(0x0100, 0x0000, 0x0000, 0x0000);
    hv507.loadRegisterToOutput();
    voltage_on_electrode1 = true;
    enabled = true;
    digitalWrite(BUTTON_OUT_PIN, HIGH);
    timerBipolarSquareWave.begin(hv507InvertOutputs, DELAY); // update time (us)
  } else {
    digitalWrite(LED_PIN, LOW);

    // Comment out for frequencies other than 0.1 Hz (so I'm logging voltages from the proper electrode)
    // timerBipolarSquareWave.end();
    // hv507.loadArrayToRegister(0x0000, 0x0000, 0x0000, 0x0000);
    // hv507.loadRegisterToOutput();
    // digitalWrite(BUTTON_OUT_PIN, LOW);
    // voltage_on_electrode1 = false;

    enabled = false;
  }
}

void hv507InvertOutputs() {
  if (voltage_on_electrode1 && enabled) {
    hv507.loadArrayToRegister(0x0000, 0x0000, 0x0000, 0x0080); // corresponds to pin 56 (56-48=8 --> plug 1000,0000 into https://www.rapidtables.com/convert/number/hex-to-binary.html?x=80)
    hv507.loadRegisterToOutput();
    voltage_on_electrode1 = false;
  } else if (!voltage_on_electrode1) {
    hv507.loadArrayToRegister(0x0100, 0x0000, 0x0000, 0x0000); // corresponds to pin 9
    hv507.loadRegisterToOutput();
    voltage_on_electrode1 = true;
  } else { // !voltage_on_electrode1 && !enabled (interrupt triggered to end voltage)
    hv507.loadArrayToRegister(0x0000, 0x0000, 0x0000, 0x0000); 
    hv507.loadRegisterToOutput();
    digitalWrite(BUTTON_OUT_PIN, LOW);
    voltage_on_electrode1 = false;
    timerBipolarSquareWave.end();
  }
}
