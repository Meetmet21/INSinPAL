INSinPAL User Guide
===================

## Table of content

[//]:

* [Introduction](#Introduction)
* [Palindromic fragile sites](#Palindromic-fragile-sites)
  * [Palindrome mining algorithm](#Palindrome-mining-algorithm)
  * [Filtering step](#Filtering-step)
  * [Interesting sites](#Interesting-sites)
* [SV callers](#SV-callers)
  * [INSurVeyor](#INSurVeyor)
  * [Manta](#Manta)
  * [Basil](#Basil)
  * [Benchmarking](#Benchmarking)
* [Annotations](#Annotations)
  * [MEIs annotation](#MEIs-annotation)
  * [Source of inserted sequence](#Source-of-inserted-sequence)
  * [Length of inserted sequence](#Length-of-inserted-sequence)
  * [Functionnal annotations](#Functionnal-annotations)
* [Installation](#Installation)  
* [Execution](#Execution)

[//]:


## Introduction

When working with Whole Genome Sequencing (WGS) data for diagnosis purposes, one has to put in place filtering steps to retrain only relevant mutations/calls made by variant callers otherwise thousands of them have to filtered by experts. INSinPAL is a Snakemake workflow analyzing, filtering and formatting large insertions calls in WGS data. INSinPAL uses two structural variant (SV) callers ([Manta](https://github.com/Illumina/manta) and [Basil](https://github.com/seqan/anise_basil)) and one insertion specialized caller ([INSurVeyor](https://github.com/kensung-lab/INSurVeyor)) to form a metacaller for germline large insertion variants detection. The main idea behind INSinPAL is to focus variant selection to palindromic fragile sites as they are prone to form DNA secondary strucures and undergo chromosomal rearrangment. Two modules allow the detection of source and length of the inserted sequence in a given breakpoint if possible which gives more information for variant analysis by experts.  


## Palindromic fragile sites

### Palindrome mining algorithm

### Filtering step

### Interesting sites


## SV callers

### INSurVeyor

### Manta

### Basil

### Benchmarking

## Annotations

### MEIs annotation

### Source of inserted sequence

### Length of inserted sequence

### Functionnal annotations


## Installation


## Execution
