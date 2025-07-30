function o = fchoose(ai,ai_1,amax)
tau1=9;
if amax<=2*ai-ai_1 
    o=amax;
else
    o=(2*ai-ai_1+min(amax,ai+tau1*(ai-ai_1)))/2;
end
end

