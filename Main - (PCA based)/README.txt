# Sinteze de HRTF a partir de dados antropométricos 
Inicialmente o método utilizado para a sinteze de HRTFs é baseado em 'shallow neural networks', que consiste em um modelo relativamente simples de arquitetura, mas que é capaz de aproximar qualquer função continua como apresentado pelo teorema universal da aproximação. Para o uso desse tipo de arquitetura o processamento dos parametros de entrada e saida são vitais, por isso aqui utilizamos da análise de componentes principais para melhorar a performance do modelo.


	## Rotinas principais
*_Rotinas ordenadas por sequencia necessária de uso._


### Preprocess_Datasets.m 
Realiza o calculo da DTF a partir dos arquivos SOFA para os dataset CIPIC, ARI, ITA, 3D3A e RIEC independentemente ou juntos. Saida DTF no domínio da frequência. Pode ser utilizado como alternativa à parte inicial de Preprocess_CIPIC.m, uma vez que a PCA ainda é necessária.  


### Anthropometry_Datasets.m
Permite escolha dos parametros antropometricos utilizados no treinamento das redes com qualquer um dos datasets, seja independentemente ou em determinados grupos. 


### PCA_PTF 
Análise de componentes principais sobre dados obtidos em Preprocess_Datasets.m. Como saída temos matrizes para treinamento da network e posterior reconstrução.

### NeuralNet_Datasets 
Criação de arquitetura e trinamento de network para cada posição de fonte relativa a cada orelha com o dataset especificado. Como saída temos os modelos treinados. 

### Rebuild_HRTF.m 
Exemplo da sintetização de HRTFs individualizadas para um indivíduo com dados desconhecidos pela network. Permite a vizualização comparativa entre hrtfs e hrirs para determinada direção entre o modelo estipulado pela network e os valores reais. Como saída temos objetos SOFA com as hrtfs simuladas e as hrtfs medidas


	## Rotinas principais (CIPIC dataset)
*_Rotinas ordenadas por sequencia necessária de uso._
*_Rotinas desenvolvidas no inicio do projeto, podem ser substituidas facilmente pelas rotinas apresentadas acima._

### Preprocess_CIPIC.m 
Realiza a tranformação: HRIR -> HRTF -> DTF para arquivos no formato CIPIC original, posteriormente a análise de componentes principais sobre os dados e preparação para target no algoritmo de treinamento da network. Além de separar os dados antropométricos correspondentes os paramentros de entrada na network. 


### NeuralNet_CIPIC.m 
Criação de arquitetura e trinamento de network para cada posição de fonte relativa a cada orelha com o dataset CIPIC. Como saída temos os modelos treinados. 



