setGlobalx(4);
testUseGlobalVar();
N = zeros(4,4);
testProc;


function setGlobalx(val)
global M;

M = val;
end

function testUseGlobalVar()
global M;
fprintf('M = %d\n', M);
end

function testProc
fprintf('M = %d\n', N(1,1));
end