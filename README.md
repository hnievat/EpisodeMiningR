# EpisodeMiningR
WINEPI Episode Mining algorithm for R

Original algorithm by Mannila, et al. 1997

Based on the Python implementation by @duytri and @nkmrtty

Uses libraries gtools and rapportools

Tested on RStudio 1.4.1717, R 4.1.1, gtools 3.9.2, and rapportools 1.0

## Features
- Program flow and variable names as close to the original as possible, with few changes
- Features both serial and parallel episode making

## Files
- WINEPI.R: the program itself
- example.csv: the example dataset in @duytri's version

## Issues
- Few sanity checks, must make sure data is in the right format

## To do
- Implement WINEPI
- Optimize code and data structure

## Thanks

- [EpisodeMining](https://github.com/duytri/EpisodeMining) by @duytri
- [episode_mining](https://github.com/nkmrtty/episode_mining) by @nkmrtty
