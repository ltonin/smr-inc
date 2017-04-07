function [psd, events, labels, settings] = smrinc_concatenate_data(dir, subject, pattern, extension)

    % Get data
    cpattern = ['*' subject '*' pattern '*'];
    Files = util_getfile(dir, extension, cpattern);
    nfiles = length(Files);
    
    psd = [];
    TYP = [];
    DUR = [];
    POS = [];
    Mk  = [];
    settings = cell(nfiles, 1);
    
    for fId = 1:nfiles
       cdata = load(Files{fId});
       info = util_getfile_info(Files{fId});
       
       % Get modality
       switch(info.modality)
           case 'offline'
               cmodality = 0;
           case 'online'
               cmodality = 1;
           otherwise
               error('chk:mod', 'Unknown modality');
       end
       
       % Concatenate events
       TYP = cat(1, TYP, cdata.events.TYP);
       DUR = cat(1, DUR, cdata.events.DUR);
       POS = cat(1, POS, cdata.events.POS + size(psd, 1));
       Mk  = cat(1, Mk, cmodality*ones(size(cdata.psd, 1), 1));
       
       % Concatenate data
       psd = cat(1, psd, cdata.psd);
        
       % Concatenate settings
       settings{fId} = cdata.settings;
    end
    
    events.TYP = TYP;
    events.DUR = DUR;
    events.POS = POS;
    labels.Mk  = Mk;
end