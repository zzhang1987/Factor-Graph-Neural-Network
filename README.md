# Factor Graph Neural Network

Created by [Zhen Zhang](https://zzhang.org), [Fan Wu](https://github.com/kkkkahlua) and [Wee Sun Lee](https://www.comp.nus.edu.sg/~leews/). 

## Requirements
The following packages are required: 

1. Python 3 
2. PyTorch 1.0
3. [AD3](https://github.com/andre-martins/AD3)

## Introduction
Most of the successful deep neural network architectures are structured, often consisting of elements like convolutional neural networks and gated recurrent neural networks. Recently, graph neural networks have been successfully applied to graph structured data such as point cloud and molecular data. These networks often only consider pairwise dependencies, as they operate on a graph structure. We generalize the graph neural network into a factor graph neural network (FGNN) in order to capture higher order dependencies. The FGNN is defined using two types of modules, the Variable-to-Factor (VF) module and the Factor-to-Variable (FV) module, as shown in Figure~\ref{fig:FGNN}. These modules are combined into a layer, and the layers can be stacked together into an algorithm. We show that the FGNN is able to exactly parameterize the Max-Product Belief Propagation algorithm, which is widely used in finding approximate \map (MAP) assignment of a PGM. Thus, for situations where belief propagation gives best solutions, the FGNN can mimic the belief propagation procedure. This repo provides the code for testing FGNN on synthetic MAP inference problem and point cloud segmentation on real dataset.

![FGNN](images/FGNN.svg?sanitize=true "Factor Graph Neural Network")

## Citation


## MAP Inference 


## Point Cloud Segmentation

Under construction
