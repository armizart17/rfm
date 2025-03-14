
generateScattering(1540,  0.05*c0*ones(100,100), nScatterers)

function map = generateScattering(propertyMap, stdMap, nScatterers)
[M,N] = size(propertyMap);
scatterProb = nScatterers/M/N;
scatterMap = rand(M,N)<scatterProb;
randMap = propertyMap + stdMap.*randn(M,N);
map = propertyMap;
map(scatterMap) = randMap(scatterMap);
end
