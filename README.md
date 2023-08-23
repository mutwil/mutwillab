# Welcome to Mutwil Lab!
Here we maintain tools to study plant biological processes and their evolution through the lens of systems biology.
Please visit our [website](https://www.plant.tools/) for more details 
## Setup
### Clone repository into local machine
```
$ git clone <placeholder>
```
### Create environment, install required packages

```
$ cd mutwillab
$ virtualenv -p python3 .venv
$ source ./.venv/bin/activate
$ pip install --upgrade pip
$ pip install -r ./setup/requirements.txt
```
That's it! You should have everything you need to import tools into your python scripts

## Importing tools
Use this boiler plate at the start of your scripts.

```
#Boiler plate
import sys
sys.path.insert(0, "/path/to/mutwillab/src")

#Import modules from mutwillab
from coexpression import pearson

#your code
```
**For Mutwil Lab members:** If writing executable script for your pipline in `main/` folder, your can use this boiler plate.
```
#!/usr/bin/env python
#Boiler plate
import sys
if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)
```
Happy coding!

# ATTENTION
This page is still being worked on. Please be patient with us. More will be revealed. :)
