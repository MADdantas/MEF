# MEF
This repository contains methodes to calculate the displacement and reaction forces on trusses and gantries
## Requirements
* Install python 3 on your computer. You can download it [here](https://www.python.org/downloads/).

* Install the numpy library, we will use it to generate matrices and do basic operations with them.

  - To install it you just need to open a CMD window and paste the following comand (do this after completing the python installation)
```cmd
pip install numpy
```

## How To Use

To solve a new problem, just create a python script in the same directory and import the function you wanna use.

For trusses:
```python
from trelica import trelica
```

For gantries:
```python
from portico import portico
```

After that import the numpy library as follows:
```python
import numpy as np
```
Then you can use the function and it will return three matrices.

The first one is the nodal displacement, the second is the nodal forces, and the last one is the nodal forces in each beam.

For example:
```python
from trelica import trelica
import numpy as np

[U,F,f] = trelica(EA,L,theta,n,b,Fe)
print("Deslocamentos nodais global: \n",U,"\n","Forças nodais: \n",F,"\n","Forças nodais em cada barra: \n",f)
```

Given the inputs, the code will output U, F, and f (The first, sencond and tird matrices as explain before).

You can find examples of how to fill the inputs in the files:
  - main_trelica.py (trusse)
  
  ![image](https://user-images.githubusercontent.com/24358380/190239151-b1f4ed01-4aea-481b-a507-8300b7854bce.png)

  - main_portico.py (gantry)
  
  ![image](https://user-images.githubusercontent.com/24358380/190239232-50ae2bcc-5994-4fb1-81da-729d04f5bb6f.png)
