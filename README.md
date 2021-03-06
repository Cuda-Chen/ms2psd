# ms2psd
Power Spectral Density Probability Density Functions Calculation
describe by [McMarana 2004](https://pubs.usgs.gov/of/2005/1438/).

## How to Install
### Requirements
- libmseed
- FFTW

### Compile
```
$ make
```

## Usage
```
./ms2psd [f1] [f2] [f3] [f4] [totype] [input] [resp]

Input parameters:
f1, f2, f3, f4: four-corner frequencies (Hz)
totype: specify the following numbers for output waveform format:
        0: displacement
        1: velocity
        2: acceleration
input: input waveform. Should be miniSEED format
resp: response file in SACPZ format
```

## TO-DO
- [x] first workable version
    - [x] PDF calculation
    - [x] PDF plotting feature
        - [x] GMT w/ bash script
        - [x] Matplotlib
- [ ] CUDA support
- [ ] PDF file description
- [ ] PDF file search utility for PDF aggregation
