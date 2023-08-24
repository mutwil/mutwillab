# Welcome to Mutwil Lab!
Here we maintain tools to study plant biological processes and their evolution through the lens of systems biology.

Please visit our [website](https://www.plant.tools/) for more details 
## Setup
### Clone repository into local machine
```
$ git clone https://github.com/mutwil/mutwillab.git
```
### Create environment, install required packages

```
$ cd mutwillab
$ virtualenv -p python .venv
$ source ./.venv/bin/activate
$ pip install --upgrade pip
$ pip install -r ./setup/requirements.txt
```
if the install failed then probably the python version used by virtualenv is not up to date. A quick fix is to install conda and use it to create an environment with python>3.10.

```
$ conda create --name python_3_10 python=3.9
$ conda activate python_3_10
$ cd mutwillab
$ virtualenv -p python .venv
$ conda deactivate
$ source ./.venv/bin/activate
$ pip install --upgrade pip
$ pip install -r ./setup/requirements.txt
```

That's it! You should have everything you need to import tools into your Python scripts


## Importing modules
If using IDEs (like VScode), use this boiler plate at the start of your scripts.

```
#Boiler plate
import sys
sys.path.insert(0, "/path/to/mutwillab/src/")
sys.path.insert(0, "/path/to/mutwillab/.venv/lib/python3.11/site-packages/")

#Import modules from mutwillab
from coexpression import pearson

#your code
```
**For Mutwil Lab members:** If writing executable script for your pipline in `main/` folder, your can use this boiler plate.
```
#Boiler plate
import sys
if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         sys.path.insert(0, parent_module)
#Import modules from mutwillab
from coexpression import pearson
```
Activate virtual enviroment and run the script from CLI.
```
$ source ./venv/bin/activate
$ python ./src/main/yourscript.py
```
Happy coding!

# ATTENTION
This page is still being worked on. Please be patient with us. More will be revealed. :)
