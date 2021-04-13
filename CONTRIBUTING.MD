Style guide

flake8 is used to check the python code against the style conventions of PEP 8:
flake8 my_notebook_name.py


Also, we use black that is a python code auto-formatter:
black -l 79 -t py36 my_notebook_name.py


We use isort to sort import statements:
isort my_notebook_name.py

Code quality
We use pylint to further check bugs and quality in the python script:
pylint my_notebook_name.py

The output of pylint may serve as input for review discussions, as seen fit.