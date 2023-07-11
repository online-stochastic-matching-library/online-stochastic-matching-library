# An Open Library of Online Stochastic Matching Algorithms

This repo contains the source code of Online Stochastic Matching algorithms related to the paper
"**Online Matching Meets Sampling Without Replacement**"

Supported Settings:
- Online Bipartite Matching 
- Online Stochastic Matching 

Supported Algorithms:
- Min Degree
- Ranking <sup>[[9]](#9)</sup>
- Feldman et al. <sup>[[3]](#3)</sup>
- Bahmani Kapralov <sup>[[1]](#1)</sup>
- Manshadi et al. <sup>[[10]](#10)</sup>
- Jaillet and Lu <sup>[[8]](#8)</sup>
- Brubach et al. <sup>[[2]](#2)</sup>
- Balance OCS <sup>[[4]](#4)</sup>
- Correlated Sampling <sup>[[6]](#6)</sup>
- Poisson OCS <sup>[[7]](#7)</sup>
- Heaupler et al. <sup>[[5]](#5)</sup>
- Top Half Sampling <sup>[[7]](#7)</sup>
- **Stochastic SWOR**
- **Regularized Greedy** 
- **Balance SWOR**

## Repository Overview
There are several directories and files in this repo:

- `\algorithms`: contains the source code for the supported algorithms mentioned above.
- `\real_world`</font>: contains raw data of large real-life graphs from Network Data Repository <sup>[[11]](#11)</sup>.
- `\real_world_result`: contains experimental results of algorithms on bipartite graphs induced from large real-life instance.
- `\real_world_small`: contains raw data of small real-life graphs from Network Data Repository <sup>[[11]](#11)</sup>.
- `\real_world_small_result`: contains experimental results of algorithms on bipartite graphs induced from small real-life instance.
- `stochastic_matching_library.cpp`: the main executive cpp file.
- `read_file.cpp`: read raw graph data from real-life instances and transfer into bipartite versions .
- `natural_lp.h`: the code implementation of ellipsoid method for solving natural LP <sup>[[6]](#6)</sup>.
- `graph.h`: the code implementation of type graph and realization graph.
- `decomposite_graph.h`: the code implementation of graph decomposition into blue and red edges in Feldman et al.<sup>[[3]](#3)</sup> and Bahmani and Kapralov <sup>[[1]](#1)</sup>.
- `cycle_break_graph.h`: the code implementation of cycle break for type graph in Jaillet and Lu <sup>[[8]](#8)</sup> and Brubach et al. <sup>[[2]](#2)</sup>.
- `flow_graph.h`: the code implementation of maximum flow by Dinic's algorithm <sup>[[12]](#12)</sup>.


## Compiling in Windows

Our library requires GNU Linear Programming Kit (GLPK) to solve linear programming.
On Windows System, GLPK for MinGW-w64 available at https://packages.msys2.org/base/mingw-w64-glp.

You may compile this library by: 
```
g++.exe stochastic_matching_library.cpp -static -O2 -lglpk -std=c++11 -Wall -o stochastic_matching_library.exe
```


## References
<a id="1">[1]</a> 
Bahman Bahmani and Michael Kapralov. Improved bounds for online stochastic matching. In
Proceedings of the 18th Annual European Symposium on Algorithms, pages 170–181. Springer,
2010. 

<a id="2">[2]</a> 
Brian Brubach, Karthik Abinav Sankararaman, Aravind Srinivasan, and Pan Xu. New algo-
rithms, better bounds, and a novel model for online stochastic matching. In Proceedings of
the 24th Annual European Symposium on Algorithms. Schloss Dagstuhl-Leibniz-Zentrum fuer
Informatik, 2016. 

<a id="3">[3]</a> 
Jon Feldman, Aranyak Mehta, Vahab Mirrokni, and Shan Muthukrishnan. Online stochastic
matching: beating 1 − 1/e . In Proceedings of the 50th Annual IEEE Symposium on Foundations
of Computer Science, pages 117–126. IEEE, 2009. 

<a id="4">[4]</a> 
Ruiquan Gao, Zhongtian He, Zhiyi Huang, Zipei Nie, Bijun Yuan, and Yan Zhong. Improved
online correlated selection. In Proceedings of the 62nd Annual IEEE Symposium on Founda-
tions of Computer Science, pages 1265–1276. IEEE, 2022. 

<a id="5">[5]</a> 
Bernhard Haeupler, Vahab S Mirrokni, and Morteza Zadimoghaddam. Online stochastic
weighted matching: Improved approximation algorithms. In Proceedings of the 7th Inter-
national Workshop on Internet and Network Economics, pages 170–181. Springer, 2011. 

<a id="6">[6]</a> 
Zhiyi Huang and Xinkai Shu. Online stochastic matching, Poisson arrivals, and the natural
linear program. In Proceedings of the 53rd Annual ACM SIGACT Symposium on Theory of
Computing, pages 682–693, 2021. 

<a id="7">[7]</a> 
Zhiyi Huang, Xinkai Shu, and Shuyi Yan. The power of multiple choices in online stochas-
tic matching. In Proceedings of the 54th Annual ACM SIGACT Symposium on Theory of
Computing, pages 91–103, 2022. 

<a id="8">[8]</a> 
Patrick Jaillet and Xin Lu. Online stochastic matching: New algorithms with better bounds.
Mathematics of Operations Research, 39(3):624–646, 2014. 

<a id="9">[9]</a> 
Richard M Karp, Umesh V Vazirani, and Vijay V Vazirani. An optimal algorithm for on-
line bipartite matching. In Proceedings of the 22nd Annual ACM Symposium on Theory of
Computing, pages 352–358, 1990. 

<a id="10">[10]</a> 
Vahideh H Manshadi, Shayan Oveis Gharan, and Amin Saberi. Online stochastic matching:
online actions based on offline statistics. Mathematics of Operations Research, 37(4):559–573,
2012. 

<a id="11">[11]</a> 
Ryan A. Rossi and Nesreen K. Ahmed. The network data repository with interactive graph
analytics and visualization. In AAAI, 2015 URL https://networkrepository.com. 

<a id="12">[12]</a> 
E.A. Dinic, Algorithm for solution of a problem of maximum flow in a network with power estimation, Soviet
Math. Doll. 11 (5), 1277-1280, (1970). (English translation by RF. Rinehart).  
