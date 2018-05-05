
Rodar frequencia a partir de um arquivo xyz

./zerg.exe frequency [file]


LISTA DE TODAS AS OPCOES DO GaInput.txt

%                                      - comentario
seed = n                               - semente
restart = yes/no                       - recomecar do arquivo restart.txt
numbert_of_cores = n                   - numero de processadores no input do gamess
printing_debug                         - 0 nada ; 1 histogram ; 2 all
project_name = text                    - nome q sera usado nos inputs do gamess
population_size = n                    - populacao, precisa ser multipla de 4
maxinum_number_of_generations = n      - numero maximo de geracoes do GA
number_of_predated = n                 - numero de individuos para serem eliminados por geracao
highlander_initial_energy = R          - numero considerado infinto
highlander_max_iteration = n           - numero de geracoes q o highlander precisa sobreviver para considerar o calculo terminado
highlander_threshold = R               - limite de tolerancia para o highlander ser deposto.
large_energy_value = R                 - parametro do auto-adjust
maximum_creation_operators_boost = R   - parametro do auto-adjust
predator_method = n                    - metodo de predacao
mutation_value = R                     - valor maximo do passo da mutacao
crossover_weight = R                   - parametro do crossover
crossover_probability = R              - parametro do crossover
sCartesianDisplacementOperators = R      - parametro do operador cartesian
alfaMinGeometricDisplacement = R       - parametro do GCDO
alfaMaxGeometricDisplacement = R       - parametro do GCDO
wGeometricDisplacement = R             - parametro do GCDO
tetaMinTwistOperator = R               - parametro do Twist
tetaMaxTwistOperator = R               - parametro do Twist
contractionMinMoveToCenter = R         - parametro do move to center
contractionMaxMoveToCenter = R         - parametro do move to center
%
%
number_of_atoms = n                    - numero de atomos
% apontar os tipos de atomos
  logo depois. Default todos-1
n_atom_type_1 = n                      - numero de atomos do tipo 1
n_atom_type_2 = n                      - numero de atomos do tipo 2
n_atom_type_3 = n                      - numero de atomos do tipo 3
n_label_type_1 = text                  - label do atomo 1
n_label_type_2 = text                  - label do atomo 2
n_label_type_3 = text                  - label do atomo 3
%
number_of_parameters = n               - ATENCAO - a lista de parametros precisa estar em seguida
A1                                     - formato - 1-1, 1-2, 1-3, 2-2, 2-3, 3-3
eps1                                   - se nao existir o tipo correspondente, elimine os parametros dele da lista
p1                                     - gupta - ordem:  A, zeta, p, q e r0.
q1
r01
A12
eps12
...
12
...
A22
eps22
...
%
%
%
user_defined_method = n                - operadores fixos - consultar codigo
gamma_creation_radius = R              - parametro da IMIGRACAO
activateIntoBfgs = n                   - make comparisson between bfgs steps
similarityMethod = n                   - (-1) no similarity; (0) distance; (1) marques rmsd
similarityDebugLevel = n               - debug level - 0(none) to 2(full)
tolSimilarity = R                      - similarity criterion cut
energyReturnBfgs = -1.0e99             - parameter for debug considerations
removeSimilarStructures = n            - 0 nao remove ; 1 remove estruturas
bestIndividualSize = n                 - numero de isomeros de comparacao
radius_factor = R                      - parametro da IMIGRACAO
max_distance_between_atoms = R         - impede distancia maxima
iterations_to_repeat_if_it_is_similar  - similaridade
interaction_potential = text           - potencial de interacao opcoes: lennardjones ; gupta ; mopac ; gamess ; MopacGamess
generation_change_interaction = n      - da geracao n para frente o potencial de interacao sera o gamess
gamess_executable_path = text          - caminho do gamess
gamess_scr_path = text                 - caminho do scr 
%
%
%
bases_files_number = n                 - numero de bases
[ ATENCAO - OS NOMES DA BASES DEVEM VIR LOGO EM SEGUIDA
base = text
base = text
...
(elas devem estar na pasta 'auxFiles' - nao e necessario incluir o caminho da pasta)
]
%
%
%
gamess_header_file = text              - gamess header (lido sempre na pasta auxFiles)


FORMATO DA BASE

-> X                - simbolo do atomo
->                  - base 
->                  - base
->    ...           - base
 $END ou $ECP       - se $ECP em baixo tem que ter o ecp se nao esse e desativado
->                  - ecp
->                  - ecp
$END

COMANDOS UTEIS DO GIT

primeira instalacao:
1.
git clone https://github.com/zerg-ga/zerg.git
2.
make


atualizando a versao
1.
git pull origin master
2.
make clean
3.
make

retorna todos os arquivos para a versao original
git reset --hard


PARA QUE O BFGS-SIMILARITY FUNCIONE É NECESSARIO INCLUIR AS SEGUINTES
LINHAS NO ARQUIVO optimize.h DA BIBLIOTECA DLIB

declaracao antes do loop (perto da linha 200)
		T xBfgs;

perto da linha 227 do arquivo optimize.h
			double dummy = f(x);
			double bfgsTest = f(xBfgs);
			if (bfgsTest < min_f)
				f_value = bfgsTest;





ANOTACOES
/*

alteracoes do qga
- portabilildade no fitness.cpp e no makefile (isnan)

proximo paper
- o highlander sai passando pela populacao com um raio maior e subindo energias - algumas vezes eliminar não os piores, mas o highlander retirar os proximos a ele. o tempo e depois que o highlander comecar a se repetir muito.

- fazer um estudo sobre o hidrogenio, avaliar os efeitos da base e etc, escolher uma base que presta para descrever os clusters.

estudo das formas de gerar clusters iniciais:
- Tem o de Sao Carlos, esfera, cubo e árvore.
- usar a algebra de polya para contar isomeros lennard jones dos bimetalicos,
  gerar todos e avaliar suas diferencas.
- colocar const no cluster operators por seguranca

- pensar sobre o highlander max iteration está 1.0e99. (CHECAR DBL_MAX)

- pop_size tem que ser multiplo de 4.

- Caso o Marques enantiomers nao convirja o rmsd minimo e usado e ele nao faz nenhum aviso para isso.

- a energia minima possibel no fitness e: -1.0e90. No Similarity::checkIntoBfgs, portanto, ele retorna -1.0e99 para parar.

- se o individuo possuir rate inferior a 0.001 ele nao e considerado para ser criado na etapa de criacao.

*/





