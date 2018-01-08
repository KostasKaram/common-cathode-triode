# common-cathode-triode

Python implementation of a simple audio tube amplifier stage.

## About

This is a digital model of a simple tube preamplifier. The implementation follows the Discrete K-method, as presented in the academic work of David Te Mao Yeh. The tube model is ripped from a DAFX paper by Kristjan Dempwolf and Udo Zoelzer.
The specific schematic modelled here is the one found on page 130 of Yeh's thesis, which includes the parasitic capacitance. Please check the **References** section for further information.

## How to use

  * You can edit the **default_params.json** file to change circuit values, processing parameters and file paths
  * Run **cct_plot.py** to show a plot of output voltages for a sine input (frequency & duration is specified in the settings json)

## TODO

 * Add a high-pass filter to block DC for proper audio use 

## References

* [Digital Implementation of Musical Distortion Circuits by Analysis and Simulation](https://ccrma.stanford.edu/~dtyeh/papers/DavidYehThesissinglesided.pdf)
* [Automated Physical Modeling of Nonlinear Audio Circuits for Real-Time Audio Effects -- Part II: BJT and Vacuum Tube Examples](https://ccrma.stanford.edu/~dtyeh/papers/yeh12_taslp.pdf)
* [A Physically-motivated Triode Model for Circuit Simulations](http://recherche.ircam.fr/pub/dafx11/Papers/76_e.pdf)

