function ptno=mrn_to_ptno(mrn)
% Input MRN string, Output Patient # string

switch str2num(mrn)
    case 957098 %1st Patient's MRN
        ptno='01'; % The 1st patient
    case 1 
        ptno='957098'; % finding the reverse

    case 1000086
        ptno='02';
    case 2
        ptno='1000086';

    case 1007969
        ptno='03';
    case 3 
        ptno='1007969';

    case 1008684
        ptno='04';
    case 4
        ptno='1008684';

    case 1024274
        ptno='05';
    case 5 
        ptno='1024274';

    case 1028588
        ptno='06';
    case 6 
        ptno='1028588';

    
    case 1029145
        ptno='07';
    case 7 
        ptno='1029145';

   
    case 1031525
        ptno='08';
    case 8 
        ptno='1031525';

    case 1035444 % Heart trouble
        ptno='09';    
    case 9 
        ptno='1035444';


    case 1044389
        ptno='10';   
    case 10 
        ptno='1044389';


    case 1055892 % Sylvian fissure complication
        ptno='11';   
    case 11 
        ptno='1055892';


    case 1069507 % Head mvmt
        ptno='12';  
    case 12 
        ptno='1069507';


    case 1072129
        ptno='13';  
    case 13 
        ptno='1072129';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% For the AQ, NNL, Olea, QT comparison - borrowed from earlier study
    case 847879
        ptno='anoq10';  
    case 400474
        ptno='anoq11';
    case 772248
        ptno='anoq12';    
    case 885665
        ptno='anoq13';   
    case 891609
        ptno='anoq14';
    case 745310
        ptno='anoq15';    
    case 891477
        ptno='anoq16'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 1072634
         ptno='14';
    case 14
        ptno='1072634';


    case 1075421
         ptno='15';  
    case 15 
        ptno='1075421';

  
    case 1085919
         ptno='16'; 
    case 16 
        ptno='1085919';


    case 1080775
         ptno='17';
    case 17 
        ptno='1080775';


%     case 
%         ptno='18';    
%     case 
%         ptno='19';   
%     case
%         ptno='20';
%     case 
%         ptno='21';    
%     case 
%         ptno='22';  
%     case 
%         ptno='23';
%     case 
%         ptno='24';    
%     case 
%         ptno='25';   
%     case 
%         ptno='26';
%     case 
%         ptno='27';    
%     case 
%         ptno='28'; 
%     case 
%         ptno='29';
%     case 
%         ptno='30';    
%     case 
%         ptno='31';   
%     case 
%         ptno='32';
%     case 
%         ptno='33';    
%     case 
%         ptno='34';  
%     case 
%         ptno='35';
%     case 
%         ptno='36';    
%     case 
%         ptno='37';   
%     case 
%         ptno='38';
%     case 
%         ptno='39';    
%     case 
%         ptno='40'; 
%     case 
%         ptno='41';
%     case 
%         ptno='42';    
%     case 
%         ptno='43';   
%     case
%         ptno='44';
%     case 
%         ptno='45';    
%     case 
%         ptno='46';  
%     case 
%         ptno='47';
%     case 
%         ptno='48';    
%     case 
%         ptno='49';   
%     case 
%         ptno='50';



% NNL, Olea, QT comparison patients
    case 847879
        ptno='01';
    case 400474
        ptno='03';
    case 772248
        ptno='05';
    case 885665
        ptno='06';
    case 891609
        ptno='07';
    case 745310
        ptno='08';
    case 891477
        ptno='09';

    otherwise
        ptno=[];
end


end
