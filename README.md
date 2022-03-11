# &#9883; Dynamique moléculaire sur une molécule triatomique &#9883;

Implémentation  d'une dynamique moléculaire sur une molécule triatomique en se basant sur le code python de l'oscillateur harmonique. On essaiera d’implémenter trois oscillateurs harmoniques couplés (mimant une molécule cyclique comme le cylopropane), ou deux liaisons et un angle (comme une molécule d’eau).  

Vous pouvez lancer le programme de dynamique moléculaire d'une molécule triatomique avec Binder.  

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/w2994a/molecular_dynamics_project/HEAD) 

---
## Utilisation

```
python triatomic_model.py file [-h] [-l] [-N NB_ITER] [-T TEMPERATURE] [-g] [-s] [-m]

positional arguments:
  file                  File name of MD result

options:
  -h, --help            show this help message and exit
  -l, --launch          Launch MD
  -N NB_ITER, --nb_iter NB_ITER
                        Number of iteration for MD (default: 10000)
  -T TEMPERATURE, --temperature TEMPERATURE
                        Temperature of model in Kelvin (default: 300)
  -g, --generate_graph  Generate graph of MD
  -s, --show            Show graph (not functional with binder)
  -m, --movie           View movie of MD
```

- L'argument positionel `file` est obligatoire : il s'agit soit du nom de sortie des résultats de la dynamique moléculaire (MD), soit des résultats de la MD à analyser.  

- Il est obligatoire d'utiliser au minimum une option parmis `-l`, `-g`, `-s` ou `-m`.  

- Les options du programme sont:
  - `-h`: affiche l'aide du programme.
  - `-l`: Lance une MD et écrit les résultats dans le fichier spécifié par l'argument positionel. Si le fichier n'éxiste pas, celui-ci est créé. (&#9888; Si le ficher existe celui-ci est écrassé).
  - `-N`: Permet de spécifier le nombre d'itération de la MD (par défaut elle est de 10000, ce qui représente 5 ps de dynamique).
  - `-T`: Permet de spécifier la température du système en Kelvin (par défaut elle est de 300 Kelvin).  
  - `-g`: Permet le génération et la sauvegarde des graphiques suivants (&#9888; les graphiques générés sont réalisés pour 1000 itérations de MD):
    - les energies potentielles, cinétiques et totales en fonction du temps de la MD (`energy.png`).
    - LA température (en Kelvin) en fonction du temps de la MD (`temperature.png`).
    - les forces appliquées sur les atomes en fonctions du temps de la MD (`force.png`).
    - Les accélérations appliquées sur les atomes en fonctions du temps de la MD (`acceleration.png`).
    - La vélocité appliquée sur les atomes en fonctions du temps de la MD (`velocity.png`).
    - Les positions des atomes en fonctions du temps de la MD, ainsi que les positions du barycentre de la molécule (`position.png`). cad. la position x de l'atome B et la position x,y de l'atom C [la postion y de l'atome B = 0 et les positions x,y de l'atome A = (0, 0)].
    - L'aire de la molécule en fonction du temps de la MD (`area.png`).
  - `-s`: Permet à l'utilisateur de générer les graphiques des variables qu'il souhaite.
  - `-m` : Permet de visualiser la dynamique moléculaire du système triatomique.

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