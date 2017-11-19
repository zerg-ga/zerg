Rodar frequencia a partir de um arquivo xyz

./zerg.exe frequency [file]




- void Creation::initialize_creation(int pop_size,
  -  Gera os operadores inicias com creation_rate (1/numero de operadores).

- void Creation::set_creation_methods(Predator &pred)
  - a escolha dos metodos e aleatoria, entao se algum n der a porcentagem
    nao e necessariamente o ultimo
    - sao 10 -> 10%
    - pop = 32
    - 10% de 32 = 3,2. Ele cria 4 para aquele operador.
    - como são 10 - apenas 8 operadores vao criar individuos.
    - mas a escolha desses 8 é aleatória.


AutoAdjust [max creation variation] [energy range]

