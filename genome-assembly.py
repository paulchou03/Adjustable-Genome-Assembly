import sys 
import pandas as pd
from typing import List, Dict, Iterable
from collections import defaultdict, deque

k = 12
max_mismatches = 3
low_thresh = 4
high_thresh = 49
step = 25

readList = []
with open('project2b_reads.fasta', 'r') as file:
    lines = file.readlines()
    for line in lines:
        line = line.strip()
        if not line.startswith(">"):
            readList.append(line)

kmerMap = {}
for read in readList:
    for i in range(len(read) - k + 1): 
        kmer = read[i:i+k]
        if kmer not in kmerMap:
            kmerMap[kmer] = 1
        else:
            kmerMap[kmer] += 1

ans = []

def isThereMismatch(genomesubstr: str, readseg: str) -> bool:
    if len(genomesubstr) != len(readseg):  
        return False 
    counter = 0
    for i in range(len(genomesubstr)):
        if genomesubstr[i] != readseg[i]:
            counter += 1
    if counter <= max_mismatches:
        return True
    else:
        return False

for kmer, count in kmerMap.items():
    if low_thresh <= count <= high_thresh:
        true_count = 1
    elif count > high_thresh:
        increments = (count - high_thresh) // step
        true_count = 2 + increments
    else:
        continue
    ans.append([kmer, true_count])

kmerList = []
for kmer, true_count in ans:
    for _ in range(true_count):
        kmerList.append(kmer)

kmerLength = len(kmerList[0])
kmers = sorted(kmerList)
ans = {}
for kmer in kmers:
    sub = kmer[:-1]
    if sub not in ans:
        ans[sub] = [kmer[1:]]
    else:
        ans[sub].append(kmer[1:])

in_degree = defaultdict(int)
out_degree = defaultdict(int)
allnodes = set(ans.keys()) 
for neighbors in ans.values():  
    for neighbor in neighbors:  
        allnodes.add(neighbor)  

for node in allnodes:
    out_degree[node] = 0  
    if node in ans:
        out_degree[node] = len(ans[node])  
        neighbors = ans[node]  
    else:
        neighbors = []  
    for neighbor in neighbors:
        in_degree[neighbor] += 1

start, end = None, None
for node in allnodes:
    diff = out_degree[node] - in_degree[node]
    if diff == 1:
        if start is None:
            start = node
        else:
            break
    elif diff == -1:
        if end is None:
            end = node
        else:
            break

if start is None:
    for node in allnodes:
        if out_degree[node] > 0: 
            start = node  
            break 

stack = [start]
path = deque()
localG = {}
for node in allnodes:
    if node in ans:  
        neighbors = ans[node]  
    else:
        neighbors = []  
    localG[node] = deque(neighbors)  

while stack:
    while localG[stack[-1]]:
        stack.append(localG[stack[-1]].popleft())
    path.appendleft(stack.pop())

genome = path[0]
for kmer in list(path)[1:]: 
    genome += kmer[-1]  

ansList = []
readLength = 48
counter = 0
for i in range(len(genome) - readLength + 1):
    for read in readList:
        readsegment = read[:48]
        genomesegment = genome[i:i+readLength]
        if isThereMismatch(genomesegment, readsegment):
            ansList.append([read, readList.index(read)])
            break

with open("output.txt", "w") as f:
    for group in ansList:
        read, index = group
        f.write(f">read_{index}\n{read}\n")
