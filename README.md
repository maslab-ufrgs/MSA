# Method of successive averages
Python implementation of the method of successive averages (MSA) for traffic assignment.

## Requirements
 * [Python 3](https://www.python.org/downloads/)
 * [Python Mathematical Expression Evaluator](https://pypi.python.org/pypi/py_expression_eval)
 * It's __not recommended__ use [Python 2.7](https://www.python.org/downloads/) due divergents results between Python 3 and Python 2 executions

 ## Networks
 Available at:  [Networks](https://github.com/maslab-ufrgs/transportation_networks)

## Usage

```sh
python3 sucessive_averages.py [OPTIONS]
```

All the options have usable defaults so check them before running an experiment.

Use:

```sh
python3 sucessive_averages.py -h
```

## Results
The results of each experiment are available in the folder "results" after you ran the experiment.

## Options

```
arguments:
  -h, --help            show this help message and exit
  -f FILE               The network file. (default: None)
  -i ITERATIONS, --iterations ITERATIONS  Number of iterations. (default: 1000)
```
