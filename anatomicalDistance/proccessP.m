%%

% p = pVals{1}*length(pVals{1});

% [~,p] = bonferroni_holm(pVals{1});

% [~,p] = bonferroni_holm(mean(pGlobal,1));

p = pVals{1};


sts = {'NS','*','**','***'};
cts = [0.05,0.01,0.001];

T = cell(size(p));

for ii = 1:length(p)
    for jj = length(cts):-1:1
        if p(ii)<cts(jj)
            T{ii} = sts{jj+1};
            break
        elseif jj==1
            T{ii} = sts{1};
        end
    end
end

B = sprintf(num2str(p,'%.1e, '));



