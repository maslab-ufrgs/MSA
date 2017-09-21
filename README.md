# Method of successive averages
Python implementation of the method of successive averages (MSA) for traffic assignment.

## Requirements
 * [Python 3](https://www.python.org/downloads/)
 * [Python Mathematical Expression Evaluator](https://pypi.python.org/pypi/py_expression_eval)
 
 ## Networks
 Available at:  [Networks](https://github.com/maslab-ufrgs/network-files)
 
## Usage

```sh
python3 sucessive_averages.py [OPTIONS]
```
Or:
```sh
./sucessive_averages.py [OPTIONS]
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
  -e EPISODES, --episodes EPISODES  Number of episodes. (default: 1000)
```
