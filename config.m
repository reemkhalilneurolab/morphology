function  flag = config()
        currPath = strsplit(pwd,'\');
        currDir = strcmp(currPath{end},'morphology')
        if currDir == 1
            baseFolder = pwd;
            cd (baseFolder)
            userpath (baseFolder)
            addpath(genpath(baseFolder)) 
            flag = 1;
        else       
           f = msgbox('please cd into morphology folder ','Error');
           flag =0;
        end    
end

