Rodar frequencia a partir de um arquivo xyz

./zerg.exe frequency [file]



LISTA DE TODAS AS OPCOES DO GaInput.txt

%                                      - comentario
seed = n                               - semente
restart = yes/no                       - recomecar do arquivo restart.txt
numbert_of_cores = n                   - numero de processadores no input do gamess
project_name = text                    - nome q sera usado nos inputs do gamess
population_size = n                    - populacao, precisa ser multipla de 4
maxinum_number_of_generations = n      - numero maximo de geracoes do GA
highlander_initial_energy = R          - numero considerado infinto
highlander_max_iteration = n           - numero de geracoes q o highlander precisa sobreviver para considerar o calculo terminado
large_energy_value = R                 - parametro do auto-adjust
maximum_creation_operators_boost = R   - parametro do auto-adjust
predator_method = n                    - metodo de predacao
mutation_value = R                     - valor maximo do passo da mutacao
crossover_weight = R                   - parametro do crossover
crossover_probability = R              - parametro do crossover
sCartesianDisplacemntOperator = R      - parametro do operador cartesian
alfaMinGeometricDisplacement = R       - parametro do GCDO
alfaMaxGeometricDisplacement = R       - parametro do GCDO
wGeometricDisplacement = R             - parametro do GCDO
tetaMinTwistOperator = R               - parametro do Twist
tetaMaxTwistOperator = R               - parametro do Twist
contractionMinMoveToCenter = R         - parametro do move to center
contractionMaxMoveToCenter = R         - parametro do move to center
number_of_atoms = n                    - numero de atomos
gamma_creation_radius = R              - parametro da IMIGRACAO
radius_factor = R                      - parametro da IMIGRACAO
max_distance_between_atoms = R         - impede distancia maxima
iterations_to_repeat_if_it_is_similar  - similaridade
interaction_potential = text           - potencial de interacao
gamess_executable_path = text          - caminho do gamess
gamess_scr_path = text                 - caminho do scr 
bases_files_number = n                 - numero de bases
[ ATENCAO - OS NOMES DA BASES DEVEM VIR LOGO EM SEGUIDA
base = text
base = text
...
(elas devem estar na pasta 'auxFiles' - nao e necessario incluir o caminho da pasta)
]
gamess_header_file = text              - gamess header (lido sempre na pasta auxFiles)






ANOTACOES
/*
- (TESTAR)se parar com o highlander, calcule a frequencia dele, se der bom, pare, se n�o SALVE esse cara e continue.

alteracoes do qga
- portabilildade no fitness.cpp e no makefile (isnan)
- CHECAR DBL_MAX

proximo paper
- o highlander sai passando pela populacao com um raio maior e subindo energias - algumas vezes eliminar n�o os piores, mas o highlander retirar os proximos a ele. o tempo e depois que o highlander comecar a se repetir muito.

- fazer um estudo sobre o hidrogenio, avaliar os efeitos da base e etc, escolher uma base que presta para descrever os clusters.

estudo das formas de gerar clusters iniciais:
- Tem o de Sao Carlos, esfera, cubo e �rvore.
- usar a algebra de polya para contar isomeros lennard jones dos bimetalicos,
  gerar todos e avaliar suas diferencas.
- colocar const no cluster operators por seguranca

- pensar sobre o highlander max iteration est� 1.0e99.
- pop_size has to be multiple of four.

*/


