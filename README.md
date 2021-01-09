# Recovery
Implementação de um método para dimensionar o diâmetro nominal de paraquedas semi-elipsoidais, os esforços atuantes, em especial o cálculo da força de abertura utilizando o método de Planz, e seus componentes (cabos de sustentação e caps), ainda resolve o sistema de equações diferenciais da dinâmica de voo do paraquedas em condição de vôos subsônicos. 

Plota gráficos da força de arrasto, velocidade e altura do sistema em função do tempo.

Plota gráficos do pico força de arrasto em função do tempo e da altura considerando a velocidade da carga de recuperação em queda livre após o apogeu.
Obs.: Os gráficos representam qual seria a força de arrasto caso a velocidade terminal fosse a do foguete em queda livre após um tempo ou distância do apogeu qualquer (Estudo em análise para verificar a coerência dos resultados).

Apresenta um gráfico da velocidade terminal por massa da carga de recuperação utilizando o paraquedas dimensionado, dando destaque ao um range de imprevisibilidade considerando a velocidade terminal máxima definida pelo cliente.

Apresenta também um gráfico com curvas de velocidade terminal por massa para vários valores de diâmetros, além do dimensionado no código. Estes valores possuem variação absoluta do diâmetro nominal cálculado de 30 cm.

Cálcula ainda a quantidade de tecido necessária para a fabricação do paraquedas.

Retorna, além dos valores dos componentes do paraquedas e dos esforços atuantes, as funções que foram utilizadas para montar cada um dos gráficos citados anteriormente

Obs.: O sistema de equações diferênciais do voo do paraquedas é apresentado no pdf Dinâmica_de_voo.
