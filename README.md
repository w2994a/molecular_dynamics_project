# &#9883; Dynamique moléculaire sur une molécule triatomique &#9883;

Implémentation  d'une dynamique moléculaire sur une molécule triatomique en se basant sur le code python de l'oscillateur harmonique. On essaiera d’implémenter trois oscillateurs harmoniques couplés (mimant une molécule cyclique comme le cylopropane), ou deux liaisons et un angle (comme une molécule d’eau).  

Vous pouvez lancer le programme de dynamique moléculaire d'une molécule triatomique avec Binder.  

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/w2994a/molecular_dynamics_project/HEAD) 

---
## Utilisation

```
python triatomic_model.py file [-h] [-l] [-n NB_ITER] [-g] [-s]


positional arguments:
  file                  File name of MD result

optional arguments:
  -h, --help            show this help message and exit
  -l, --launch          Launch MD
  -n NB_ITER, --nb_iter NB_ITER
                        Number of iteration for MD (default: 10000)
  -g, --generate_graph  Generate graph of MD analysis
  -s, --show            Show graph (not functional with binder)
```

---
## Installation sur votre pc local
1. Assurrez-vous d'avoir une installation miniconda ou anaconda (à l'aide de la commande : `conda --version`) et le logiciel de gestion de version Git. Si ce n'est pas le cas procéder à une installation de `conda` ([Miniconda](https://docs.conda.io/en/latest/miniconda.html) • [Anaconda](https://www.anaconda.com/products/individual)) et/ou de `Git` ([Git](https://git-scm.com/downloads)) 

2. Clonez le dépôt du programme de dynamique moléculaire :  
```
git clone https://github.com/w2994a/molecular_dynamics_project.git
```
ou via connexion ssh : 
```
git clone git@github.com:w2994a/molecular_dynamics_project.git
```  

3. Déplacez-vous dans le répertoire du dépôt du programme de dynamique moléculaire :
```
cd molecular_dynamics_project
```  

4. Créez l'environnement conda pour le programme de dynamique moléculaire :
```
conda env create -f binder/environment.yml
```  

5. Activez l'environnement pour le programme de dynamique moléculaire :
```
conda activate molecular-dynamics
```  

---  

&#169; Willam Amory, Lucas Rouaud - 2022. &#129418;