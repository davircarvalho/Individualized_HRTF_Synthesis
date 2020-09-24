function [PCWs, coeffs, eigen_value] = PCA_fun(HRTF_data, no_PC)
%% General info 
% Realiza a analise de componentes principais pelo metodo de decomposição 
% dos autovalores a partir da matriz de covariancia ('eig')
%%% INPUTS %%%
%   - HRTF_data: matriz (NxM) em que N corresponde a observações e M
%                variáveis independentes

%   - no_PC: número de componentes principais de interesse

%%
[no_of_samples,length_training_set,no_of_directions,no_of_channels] = size(HRTF_data);

coeffs = zeros(no_of_samples,no_PC,no_of_directions,no_of_channels);
PCWs = zeros(no_PC,length_training_set,no_of_directions,no_of_channels);

for m = 1:no_of_channels
    channel = m;
    for n = 1:no_of_directions
        direction = n;

        data_mtx = HRTF_data(:,:,direction,channel);
        cov_mtx = data_mtx'*data_mtx;
        %calculating sigular vectors and corresponding sigular values
        [eigen_vector, eigen_value_mtx] = eig(cov_mtx);
        [eigen_value, sorted_index] = sort(diag(eigen_value_mtx),'descend');
        temp = eigen_vector;
        for i = 1:length(sorted_index)
            eigen_vector(:,i) = temp(:,sorted_index(i,1));   
        end

        %eigen vectors
        eigen_vector = data_mtx*eigen_vector;
        for i = 1:length_training_set
            eigen_vector(:,i) = eigen_vector(:,i)/norm(eigen_vector(:,i));  %normalized eigenvectors
        end

        %The principal component variances are the eigenvalues 
        % of the covariance matrix of the input.
        eigen_value = (eigen_value/sum(eigen_value))*100;%normalized eigen values
        PCWs(:,:,n,m) = eigen_vector(:,1:no_PC)'*data_mtx;%low dimention projection
        coeffs(:,:,n,m) = eigen_vector(:,1:no_PC);   
    end
end
end