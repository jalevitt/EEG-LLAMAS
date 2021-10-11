% instantiate the library
disp('Loading the library...');
lib = lsl_loadlib();

% resolve a stream...
disp('Resolving an EEG stream...');
result = {};
while isempty(result)
    result = lsl_resolve_byprop(lib,'name','Keyboard'); 
end

disp('Opening an inlet...');
Keyinlet = lsl_inlet(result{1});
info = result{1}.channel_format()
disp('Now receiving chunked data...');

KeyPresses = cell(10000, 1);
TimeStamps = cell(10000, 1);
i=1;
6.541712626000000e-315;
Terminate = false;
while Terminate == false
    [mrks,ts] = Keyinlet.pull_chunk();
    if ~isempty(mrks)
        KeyPresses(i) = {mrks};
        TimeStamps(i) = {ts};
        i=i+1;
        mrks
        ts
        if mrks == 6.541712626000000e-315
            Terminate = true;
        end
    end
    pause(0.01)
end