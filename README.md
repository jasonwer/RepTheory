# RepTheory
Representaion theory - a mathematical libarary for constructing general linear (super) modules

Example usage:
* getBasis [2,0,0] 
[[[2,0,0],[2,0],[2]],[[2,0,0],[2,0],[1]],[[2,0,0],[2,0],[0]],[[2,0,0],[1,0],[1]],[[2,0,0],[1,0],[0]],[[2,0,0],[0,0],[0]]]

* weightFromGT [[4,2,1,0],[3,2,0],[2,1],[1]] 
[1,2,2,2]

* dimensionW [4,0,0] 
15

* dimensionW [5,3,3,1,0] 
1890

### Outputs

Basis patterns can be exported to a TeX file. Graphical representation of a module is achieved by exporting to a Python script that uses matplotlib (under construction).
