PROGRAMA QUE AUXILIA NO CONTROLE DOS PROGRAS QUÂNTICOS




"quantumOptions.inp" - arquivo que contêm as opcoes.





USANDO: WriteQuantumInput
        [Objeto que constrói um input para um programa
         quântico usando um conjunto de coordenadas]

  -> CONSTRUCTOR

  WriteQuantumInput(
  vector<string> options ->
  options[0] = type: "mopac" ou "mopac2009"
  options[1] = projectName

  "mopac" ou "mopac2009"
    options[2] = "cabecalho do mopac"
    options[3] = "nome dos arquivos externos"
    options[4] = "nome do átomo central"
    )


  para "gamess"
    options[2-n]   = varias linhas q indicam o cabecalho do gamess
                     n termina quando ele encontra o flag: "EndOfHeader"
    options[n-m]   = varias linhas que indicam os ARQUIVOS que contem
                     as bases dos atomos.
                     m termina quando o flag "EndOfBasis" é encontrado
    options[m+1]   = se tiver a flag ActivateEcp, continue lendo.
    options[m+1-l] = varias linhas que indicam os ARQUIVOS que contem os
                     ECP's.
                     l termina com a flag "EndOfEcp"

    ATENCAO - O nome atomico precisa estar igual na tabela periodica
              primeira letra maiuscula e as outras minusculas.

    ATENCAO - E necessario a pasta auxFiles com os arquivos citados no
              options.

    ATENCAO - Para o ECP os atomos de mesmo ecp precisam estar agrupados

  -> CRIACAO DE INPUT

  string createInput(
  vector<CoordXYZ> coordinates
  sao as coordenadas que deseja-se adicionar ao input.

  int index
  numero a ser adicionado no final, serve para manter
  os inputs caso faça-se um grande número de contas.
  )

  retorna uma string com o nome do input.
  observacao -> para executa-la e necessario adicionar
  a extensao correspondente ".mop" para mopac.

  -> FUNÇÕES PÚBLICAS

  changeProjectName(
  string newProjectName
  todos os futuros inputs criados com esse objeto terao o novo
  nome designado aqui.
  )

  changeMopacHeader(
  string newMopacHeader
  todos os futuros inputs criados terao o cabeçalho definido aqui
  util para calcular frequencia, por exemplo.
  )


USANDO: ReadQuantumOutput
        [objeto que lê as propriedades de output quântico]

  CONSTRUCTOR
  -> ReadQuantumOutput(
  string type
  "mopac" , "mopac2009" ou "gamess"
  )

  LEITURA DO OUTPUT
  -> readOutput(
  string inputName
  nome do output a ser lido. Pode ser a mesma saída do
  WriteQuantumInput::createInput
  )

  observação1: a leitura usa sempre a extensao ".out"

  observação2: sempre tenta ler tudo com  os criterios definidos
               em ReadQuantumOutput, o que ele achar, é isso mesmo.

  observação3: as informacoes ficam dentro do objeto

  PÚBLICAS

  -> activateDeactivateReadings(std::string activateOption, bool activate) 
     funcao que permite desativar certas leituras.
     activate - se for true, vai ler a opcao, se for false nao vai ler.
     activateOption:
         "coordinates"
         "energy"
         "ionization"
         "dipole"
         "frequency"

  vector<CoordXYZ> getCoordinates()
  double getEnergy()
  double getIonizationPotential()
  vector<double> getDipole()
  double getFirstFrequency()


IMPLEMENTAÇÂO

Funcoes que precisam ser implementadas:

WriteQuantumInput::createInput
WriteQuantumInput::setInputProperties

ReadQuantumOutput() constructor
-> como ler
ReadQuantumOutput::readCoordinates
ReadQuantumOutput::readEnergy
ReadQuantumOutput::readIonization
ReadQuantumOutput::readDipole
ReadQuantumOutput::readFrequency
-> quando ler
ReadQuantumOutput::haveToReadCoordinates
ReadQuantumOutput::haveToReadEnergy
ReadQuantumOutput::haveToReadIonization
ReadQuantumOutput::haveToReadDipole
ReadQuantumOutput::haveToReadFrequency




UTILIDADES

-> Em GamessCalculateFrequency -> tem ums funcoes que coletam a energia de arquivos xyz.


