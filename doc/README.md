
# Short Tutorial for LATTE

This is short tutorial on how to install and use `LATTE` for traveltime computation, FATT, source location, and joint location-tomography.

## Table of Contents
- [Installation](#installation)
- [Parameters](#parameters)
- [Examples](#examples)


## Installation

Installation of `LATTE` is straightforward:

```bash
git clone https://github.com/lanl/latte_traveltime.git
cd latte_traveltime
cd src
ruby install.rb clean
```

You need to install [`FLIT`](https://github.com/lanl/flit) before installing `LATTE`. 


## Geometry

To use `LATTE` for forward modeling or FWI, a source-receiver geometry file must be provided. The format of `LATTE`'s geometry file is:

### Overall geometry file

`LATTE` requires an overall geometry file (say, `./geometry/geometry.txt`) in the following form:
 
```ruby
	shot_1_geometry.txt
	shot_2_geometry.txt
	shot_3_geometry.txt
	....
```

where each row is the name of each common-shot gather's geometry file. The names of the individual geometry files can be arbitrary, e.g., 

```ruby
	shot_1_geometry.txt
	source_222_geometry.txt
	event_33_sr.txt
	...
```

However, these geometry files should be distinct (otherwise, it does not make much sense to include two identical common-shot gathers in one survey) and should present in the same directory with the overall geometry file, i.e., the path to individual geometry files is `./geometry/shot_1_geometry.txt`, ...

### Individual geometry file

An individual geometry file is to describe the source-receiver distribution as well as source parameters. The form is:

```ruby
shot id

number of point sources
point source 1 x y z t0
point source 2 x y z t0
...

number of receivers
receiver 1 x y z weight
receiver 2 x y z weight
...
```

For example:
```ruby
1								# The unique id of this shot is 1

1								# The shot contains 1 point source
100.0 500.0 10.0 0.0			# The x, y, z location of this point source is (100, 500, 10) meters,
								# and the origin time is 0 sec.

2								# The shot contains 2 receivers
200.0 400.0 10.0 1.0			# The 1st receiver locates at (200, 400, 10) meters, with a weight of 1
1200.0 1400.0 100.0 1.0			# The 2nd receiver locates at (1200, 1400, 100) meters, with a weight of 1
```

## Model


## Dimension

## Other parameters

<!-- > **`n1` (integer)** 
- **Description**: Number of grid points along axis-1 of the generated random model.
- **Default**: `128` -->