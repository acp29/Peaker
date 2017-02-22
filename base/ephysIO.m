% Function File: ephysIO
%
%  USAGE: To save a column-major XY array of electrophysiology data:
%         ephysIO (filename,array,xunit,yunit,names,notes)
%
%         To load column-major XY array of electrophysiology data:
%         [array,xdiff,xunit,yunit,names,notes,saved] = ephysIO (filename)
%         [array,xdiff,xunit,yunit,names,notes,saved] = ephysIO ({filename,channel})
%
%  Electrophysiology data input-output functionality for MATLAB. The file
%  format is automatically detected from the filename extension. The channel
%  number of the data can be specified in the first input argument by providing
%  the filename and channel number in a cell array as shown above.
%
%  The numeric data in the array argument is in column-major XY format:
%     The first column is the X dimension
%     All subsequent columns are the Y dimension
%  All Y dimension data must have the same units.
%  The X dimension can have constant or variable sampling intervals.
%
%  xunit and yunit arguments are characters defining the X (e.g. s) and Y
%  (e.g. A, V, au) dimension units respectively .
%
%  names is a cell array of strings with the title of each column of array.
%  notes is a cell array containing comments etc.
%
%  Read and write support is provided for Axon text files (.atf), Igor
%  text files (.itx) and ephysIO's matlab binary file format (.mat).
%
%  Tab-delimited text files (.txt) or comma-separated values text files
%  (.csv) containing the data array in the above format can also be read.
%  If a header line is present it will be detected automatically.
%
%  Zip (.zip) or gzip (.gz) compressed files of the supported input file
%  formats will automatically be decompressed.
%
%  Support for reading non-matlab binary files is available via the
%  following third party helper functions distributed with ephysIO:
%    readMeta.m from ACQ4 (Luke Campagnola)
%    IBWread.m from Jakub Bialek
%    abfload.m from Harald Hentschke, Forrest Collman and Ulrich Egert
%
%  Supported input file formats:
%    ACQ4 binary (hdf5) files (.ma)
%    Stimfit binary (hdf5) files (.h5)
%    ephysIO HDF5/MATLAB binary files (.mat)
%    Igor text files (.itx or .awav)
%    Axon text files (.atf)
%    Tab-delimited text files (.txt)
%    Comma-separated values text files (.csv)
%    Axon binary files 1 and 2 (.abf)
%    Igor binary wave files (.ibw or .bwav)
%
%  Supported output file formats:
%    ephysIO HDF5 (MAT v7.3) binary files (.mat)
%    Igor text files (.itx)
%    Axon text files (.atf)
%
%  ephysIO v1.3 (last updated: 27/07/2016)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/
%
%  Copyright 2010 Andrew Charles Penn
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.


function [array,xdiff,xunit,yunit,names,notes,saved] = ...
         ephysIO (arg1,array,xunit,yunit,names,notes)

  %% Evaluate input arguments (with error checking)
  % Check supported filetypes for loading
  if (nargin == 1)
    if iscell(arg1)
      filename = arg1{1};
      ch = arg1{2};
    else
      filename = arg1;
      ch = 1;
    end
    % Get file path
    pathstr = fileparts(filename);
    % Check file exists
    if ~exist(filename,'file')
      error('File not found')
    end
    if ~isempty(regexpi(filename(end-2:end),'.gz'))
      if exist('gunzip')
        gunzip(filename,pathstr);
        filename(end-2:end)='';
      else
        error('gunzip function is not available for file decompression')
      end
    end
    if ~isempty(regexpi(filename(end-3:end),'.zip'))
      if exist('unzip')
        unzip(filename,pathstr);
        filename(end-3:end)='';
      else
        error('unzip function is not available for file decompression')
      end
    end
    pat = ['(.\.mat)*(.\.txt)*(.\.csv)*(.\.itx)*(\.awav)*(.\.atf)*(.\.ibw)*(.\.abf)*'...
           '(.\.ma)*(.\.h5)*'];
    if isempty(regexpi(filename(end-4:end),pat))
      error('Unsupported filetype')
    end
  end
  % Check supported filetypes for saving
  if nargin > 1
    if iscell(arg1)
      filename = arg1{1};
    else
      filename = arg1;
    end
    pat = '(.\.mat)*(.\.itx)*(.\.awav)*(.\.atf)*';
    if isempty(regexpi(filename(end-4:end),pat))
      error('Unsupported filetype')
    end
    ncols = size(array,2);
    if ncols < 2
      error('The data must have X and Y dimensions')
    end
    if nargin < 3
      xunit = '';
    end
    if nargin < 4
      if ~isstr(xunit)
        error('xunit must be a character string')
      end
      yunit = '';
    end
    if nargin < 5
      if ~isstr(yunit)
        error('yunit must be a character string')
      end
      % Create default column names for data array
      names = cell(ncols,1);
      if ncols > 1
        if strcmp(xunit,'s')
          names{1} = 'Time';
        else
          names{1} = 'XWave';
        end
      else
        names{1}='YWave';
      end
      for i=1:ncols-1
        names{i+1}=sprintf('YWave%d',i);
      end
    end
    if nargin < 6
      if ~iscell(names)
        error('names must be a cell array')
      end
      if numel(names) ~= ncols
        error('names must match the dimensions of the data array')
      end
      notes={};
    end
    if nargin == 6
      if ~iscell(notes)
        error('notes must be a cell array')
      end
    end
  end

  %% Perform input-output operation depending on function usage
  % LOAD DATA if only filename input argument is provided
  if (nargin == 1)
    if strcmpi(filename(end-3:end),'.mat')
      [array,xdiff,xunit,yunit,names,notes,saved] = MATload (filename);
    elseif strcmpi(filename(end-3:end),'.txt')
      [array,xdiff,xunit,yunit,names,notes] = TXTread (filename,'\t');
    elseif strcmpi(filename(end-3:end),'.csv')
      [array,xdiff,xunit,yunit,names,notes] = TXTread (filename,',');
    elseif strcmpi(filename(end-3:end),'.itx') || strcmpi(filename(end-4:end),'.awav')
      [array,xdiff,xunit,yunit,names,notes] = ITXread (filename);
    elseif strcmpi(filename(end-3:end),'.atf')
      [array,xdiff,xunit,yunit,names,notes] = ATFread (filename);
    elseif strcmpi(filename(end-2:end),'.ma')
      [array,xdiff,xunit,yunit,names,notes] = MAload (filename,ch);
    elseif strcmpi(filename(end-2:end),'.h5')
      [array,xdiff,xunit,yunit,names,notes] = H5load (filename,ch);
    elseif strcmpi(filename(end-3:end),'.ibw') || strcmpi(filename(end-4:end),'.bwav')
      if exist('IBWread')
        [array,xdiff,xunit,yunit,names,notes] = IBWload (filename);
      else
        error('The required helper function IBWread cannot be found')
      end
    elseif strcmpi(filename(end-3:end),'.abf')
      if exist('abfload')
        [array,xdiff,xunit,yunit,names,notes] = ABF2load (filename,ch);
      else
         error('The required helper function abfload cannot be found')
      end
    end
    % Convert character array of column names to a cell array and transpose
    % for column-major order (if applicable)
    if ischar(names)
      names = cellstr(names)';
    end
    % Convert character array of notes to a cell array if applicable
    if ischar(notes)
      notes = cellstr(notes);
    end
    % Calculate the number of columns in the data array
    ncols = size(array,2);
    if ncols < 2
      error('The data must have X and Y dimensions')
    end
    % Remove unit prefixes and scale the data appropriately
    [array(:,1),xunit,SF] = scale_units(array(:,1),xunit);
    xdiff = xdiff*SF;
    if abs(array(1,1))>0 & strcmpi(filename(end-3:end),'.mat')
      warning('Timebase offset will be reset to zero','TimeOffset')
      array(:,1) = xdiff * [0:size(array,1)-1]';
    end
    [array(:,2:end),yunit] = scale_units(array(:,2:end),yunit);
    if ~exist('saved','var')
      saved='';
    end
  end

  % Make wave names comply with standard name rules
  for i = 1:ncols
    % Replace commas and blanks with underscore
    names{i} = regexprep(names{i},'[\s,]','_');
    % Remove all other illegal characters
    names{i} = regexprep(names{i},'[^A-Z_0-9]*','','ignorecase');
    if isempty(names{i})
      names{i}='?';
    end
    % Remove leading numeric or underscore characters
    while isempty(regexpi(names{i}(1),'[A-Z]'))
      names{i}(1)='';
      % If wave name is empty, assign an arbitrary wave name
      if (numel(names{i})==0)
        if i == 1
          if strcmp(xunit,'s')
            names{1} = 'Time';
          else
            names{1} = 'XWave';
          end
        elseif i > 1
          names{i}= sprintf('YWave%d',i-1);
        end
      end
    end
    % Truncate wave name to 31 characters
    names{i}(31:end)=[];
  end

  % SAVE DATA if both filename and data input arguments are provided
  if nargin > 1
    % Remove unit prefixes and scale the data appropriately
    [array(:,1),xunit] = scale_units(array(:,1),xunit);
    [array(:,2:end),yunit] = scale_units(array(:,2:end),yunit);
    % Set data values below the machine precision to 0
    array(abs(array)<eps) = 0;
    if strcmpi(filename(end-3:end),'.mat')
      if sum(isnan(array(:))+isinf(array(:)))>0
        ext=strcat('.',input(sprintf(['The ephysIO filetype cannot store NaN or Inf\n'...
                                     'Please select alternative file extension '...
                                     '(e.g. itx): ']),'s'));
        filename(end-3:end)=ext;
      end
    end
    if strcmpi(filename(end-3:end),'.mat')
      MATsave (filename,array,xunit,yunit,names,notes);
    elseif strcmpi(filename(end-3:end),'.itx') || strcmpi(filename(end-4:end),'.awav')
      ITXwrite (filename,array,xunit,yunit,names,notes);
    elseif strcmpi(filename(end-3:end),'.atf')
      ATFwrite (filename,array,xunit,yunit,names,notes);
    end
    clear array xdiff xunit yunit names notes
  end

%--------------------------------------------------------------------------

function [data, unit, SF] = scale_units (data, unit)

  % Scale the data so that their respective units are without prefixes
  % Establish key-value pairs
  key = [121,122,97,102,112,110,181,109,107,77,71,84,80,69,90,89];
  val = nonzeros(-24:3:24);

  % Read unit prefix
  if numel(unit) > 1
    idx = findstr(char(key),unit(1));
    SF = 10^val(idx);
  else
    % Unit has no prefix
    idx = [];
    SF = 1;
  end

  % Scale the data and correct the unit
  if ~isempty(idx)
    % Scale the data
    data = data * SF;
    % Remove the unit prefix
    unit(1) = '';
  else
    % do nothing
  end


%--------------------------------------------------------------------------

function MATsave (filename,array,xunit,yunit,names,notes)

  %% File format version 1 - no longer used
  %% Modify the X dimension data as the difference between X values to
  %% maintain the dynamic range when converting it to single precision
  %dx = diff(array(:,1));
  %if any(diff(dx) > 1.192093e-07)
  %  xdiff = 'variable';
  %else
  %  xdiff = 'constant';
  %end
  %array(2:end,1) = dx;

  %% Round double array to 7 significant figures
  %array(array==0) = NaN;
  %scale = 10.^floor(log10(abs(array)))*1e-6;
  %array = round(array./scale).*scale;
  %array(isnan(array)) = 0;

  %% Convert data array to single machine precision for more efficient storage
  %% and transpose to row-major order so that dimensions are consistent with
  %% the char array of names
  %array = single(array)';

  % Convert names and notes cell arrays to char arrays
  %names = char(names);
  %notes =  char(notes);

  %% Save variables in matlab binary file
  %try
  %  % Save in the compressed MATLAB v7 filetype
  %  save (filename,'array','scale','start','xdiff','xunit',...
  %                 'yunit','names','notes','-v7')
  %catch
  %  % Save in default MATLAB filetype
  %  save (filename,'array','scale','start','xdiff','xunit',...
  %                 'yunit','names','notes','-mat')
  %end

  %% File format version 2
  %% Much more efficient data storage
  % Transpose data array to row-major order so that dimensions are consistent with
  % the char array of names
  array = array';
  start = array(:,1);

  % Transform data array by representing it as difference values
  % - approximate first derivative
  array = diff(array,1,2);

  % Check sampling properties of the X dimension
  if any(diff(array(1,:)) > 1.192093e-07)
    % Variable sampling interval
    xdiff = 0;
    xname = '';
  else
    % Constant sampling interval
    if abs(start(1))>0
      warning('Timebase offset will be reset to zero','TimeOffset')
    end
    xdiff = array(1,1);
    array(1,:) = [];
    xname = names(1);
    names = names(2:end);
    start(1)=[];
  end

  % Scale each row of the transformed data array by a power-of-2 scaling factor
  % Power of 2 scaling is more computationally efficient
  nrows = size(array,1);
  for i=1:nrows
    maxval = max(abs(array(i,:)));
    scale(i,1) = fix(log2(32767/maxval));
    array(i,:) = array(i,:) * 2^scale(i,1);
  end

  % Change class of variables for more efficient data storage
  % Note that we store the power-of-2 exponent for the scale factor
  start = single(start);
  scale = uint8(scale);
  array = int16(round(array));
  xdiff = double(xdiff);
  xname = char(xname);
  names = char(names);
  notes = char(notes);
  xunit = char(xunit);
  yunit = char(yunit);

  % Create variable recording serial number of date and time
  saved = now;

  % Save variables in matlab binary file
  try
    % Save in the HDF5-based MATLAB v7.3 filetype
    save (filename,'array','scale','start','xdiff','xunit',...
                   'xname','yunit','names','notes','saved','-v7.3')
  catch
    % Save in default MATLAB filetype if v7.3 is not option is not available
    save (filename,'array','scale','start','xdiff','xunit',...
                   'xname','yunit','names','notes','saved','-mat')
  end

%--------------------------------------------------------------------------

function [array,xdiff,xunit,yunit,names,notes,saved] = MATload (filename)

  % Load variables
  load(filename);

  if isa(array,'single')

    %% File format version 1
    %% Included for backwards compatibility
    % Convert data array to double precision and transpose to column-
    % major order
    array = double(array)';

    % Round double array to 7 significant figures
    array(array==0) = NaN;
    scale = 10.^floor(log10(abs(array)))*1e-6;
    array = round(array./scale).*scale;
    array(isnan(array)) = 0;

    % Calculate X dimension
    if strcmp(xdiff,'constant')
      x0 = array(1,1);
      delta = array(2,1);
      array(:,1) = delta*[0:size(array,1)-1]'+x0;
    elseif strcmp(xdiff,'variable')
      array(:,1) = cumsum(array(:,1));
    end

  elseif isa(array,'int16')

    %% File format version 2
    % Convert classes of numeric variables and transpose to column-
    % major order
    array = double(array)';
    start = double(start)';
    scale = double(scale)';
    xname = cellstr(xname)';
    names = cellstr(names)';
    notes = cellstr(notes);
    if exist('saved','var')
      saved = datestr(saved,'yyyymmddTHHMMSS');
    else
      saved = '';
    end

    % Calculate power of 2 scale factor
    scale = 2.^scale;

    % Rescale transformed data array
    ncols = size(array,2);
    for i=1:ncols
       array(:,i) = array(:,i)/scale(i);
    end

    % Backtransform data array to real world values
    array = cat(1,start,array);
    array = cumsum(array,1);

    % Calculate X dimension for constant sampling interval
    if xdiff > 0
      x = xdiff * [0:size(array,1)-1]';
      array = cat(2,x,array);
      names = cat(2,xname,names);
    else
      % do nothing
    end

  end

%--------------------------------------------------------------------------

  function [array,xdiff,xunit,yunit,names,notes] = TXTread (filename,delim)

  % Read delimited text file
  try
    % Get the data
    array = dlmread(filename,delim,0,0);
    % Evaluate data array dimensions
    ncols = size(array,2);
    % Create blank names and notes variables
    names = '';
    notes = '';
  catch
    try
      % Assume there is a header line to read
      tmpstr = char(textread(filename,'%s',1,'delimiter','\n'));
      tmp = strread(tmpstr,'%s','delimiter',delim);
      % Get the data
      array = dlmread(filename,delim,1,0);
      % Evaluate data array dimensions
      ncols = size(array,2);
      % Evaluate header line
      if numel(tmp) == ncols
        % Assign header entries to names if the dimensions match
        names = char(tmp);
        notes = '';
      else
        % Otherwise assign the header line to the notes array
        names = '';
        notes = sprintf(tmpstr);
      end
    catch
      error('The text file could not be loaded')
    end
  end

  % Evaluate X dimension
  dx = diff(array(:,1));
  if any(diff(dx) > 1.192093e-07)
    xdiff = 0;
  else
    xdiff = dx(1);
  end

  % If names array is blank, assign an array of question marks to it
  if strcmp(names,'')
    names = char(zeros(ncols,1)+63);
  end

  % Create blank units variables
  xunit = '';
  yunit = '';

%--------------------------------------------------------------------------

function ITXwrite (filename,array,xunit,yunit,names,notes)

  % Evaluate the X dimension
  x0 = array(1,1);
  dx = diff(array(:,1));
  if any(diff(dx) > 1.192093e-07)
    xdiff = 0;
    k = 1;
  else
    delta = dx(1);
    xdiff = delta;
    k = 2;
  end

  % Open a file for writing
  fid = fopen(filename,'w');

  % Print IGOR keyword header
  fprintf(fid,'IGOR\n');

  % Print WAVES line
  if ~ischar(names)
    names = char(names);
  end
  ncols = size(array,2);
  fprintf(fid,'WAVES');
  for i=k:ncols
    fprintf(fid,strcat('\t',deblank(names(i,:))));
  end
  fprintf(fid,'\n');

  % Print data matrix
  fprintf(fid,'BEGIN\n');
  if k==1
    % X dimension: 64-bit double precision variables represent data to about
    % 15 significant figures (14 decimal places)
    dataformat = '\t%.15g';
  elseif k==2
    % Y dimension: 32-bit single precision variables represent data to about
    % 7 significant figures (6 decimal places)
    dataformat = '\t%.7g';
  end
  for i = k:ncols-1
    % Y dimension: 32-bit single precision variables represent data to about
    % 7 significant figures (6 decimal places)
    dataformat = strcat(dataformat,'\t%.7g');
  end
  dataformat = strcat(dataformat,'\n');
  fprintf(fid,dataformat,array(:,k:ncols)');
  fprintf(fid,'END');
  % Print footer
  for i = k:ncols
    if k==2
      fprintf(fid,'\nX SetScale/P x %g,%g,"%s", %s; ',...
              min(array(:,1)),delta,xunit,deblank(names(i,:)));
    else
      fprintf(fid,'\nX ');
    end
      fprintf(fid,'SetScale/I y %g,%g,"%s", %s',...
              min(array(:,i)),max(array(:,i)),yunit,deblank(names(i,:)));
  end

  % Print notes as comment lines
  nopth = numel(notes);
  if nopth > 0
    for i = 1:nopth
      fprintf(fid,char(['\nX // ',notes{i}]));
    end
  end

  % End file with a blank line
  fprintf(fid,'\n');

  % Close file identifier
  fclose(fid);

%--------------------------------------------------------------------------

function [waves,xdiff,xunit,yunit,names,notes] = ITXread (filename)

  % Open a file for reading
  fid = fopen(filename,'r');

  % Read IGOR keyword header (with error checking)
  tmp = fgetl(fid);
  tmp = strjust(tmp,'left');
  if regexpi(tmp,'IGOR','once') ~= 1
    error('bad Igor text file. (missing keyword IGOR)')
  end
  tmp(1:4) = [];
  if ~isempty(deblank(tmp))
    error('unknown keyword in Igor Text file')
  end

  %% INITIALIZE
  % Start counter for the blocks of data
  block = 0;
  % Start counter for waves in the current data block
  count = 0;
  % Create an empty cell array for the wave names
  names = cell(0);
  % Create an empty cell array for the wave notes
  notes = cell(0);
  % Create an empty cell array for the wave data
  waves = cell(0);
  % Create empty starting index variable
  startidx = [];
  % Create empty delta variable
  delta = [];
  % Create empty x-axis unit variable
  xunit = [];
  % Create empty y-axis unit variable
  yunit = [];

  % Read through the blocks of data
  while (true)
    tmp = fgetl(fid);

    if strcmp(tmp,'')
      % Skip blank line

    elseif all(tmp==-1)
      % End-of-file (EOF) condition breaks from while loop
      break

    elseif ~isempty(regexpi(tmp,'WAVES','once'))
      tmp = strjust(tmp,'left');
      % Step-up data block counter
      block = block + 1;

      %% READ WAVES LINE
      % Read WAVE keyword and associated flags (with error checking)
      tmp(1:5) = '';
      pat = '([,\s]*/\s*[^''N\s]*)*';
      [junk,eop] = regexpi(tmp,pat,'once');
      tmpstr = tmp(1:eop);
      tmp(1:eop) = [];

      % Check for unknown flags
      pat = '[^/ORCSDIWBUT\s]*';
      if ~isempty(regexpi(tmpstr,pat))
        error('unknown flag')
      end

      % Issue warnings/errors for unsupported flags
      if ~isempty(regexpi(tmpstr,'T','once'))
        error('Text waves are not supported')
      end
      pat = '[RCSIWBU]*';
      if ~isempty(regexpi(tmpstr,pat))
        warning(['Ignoring WAVES flags.'...
                 ' Data loaded as double precision floating point.'])
      end

      % Read wave dimensions flag
      if strcmpi(tmp(1),'N')
        pat = 'N\s*=\s*\(.+\)';
        [junk,eop] = regexpi(tmp,pat,'once');
        tmpstr = tmp(1:eop);
        tmp(1:eop) = [];
        pat = 'N\s*=\s*\(\s*\d+\s*(,\s*\d+\s*)*)';
        if isempty(regexpi(tmpstr,pat,'once'))
          error(['WAVES declaration error. '...
                 ' Improper specification of wave dimension sizes'])
        end
        pat = 'N\s*=\s*\(\s*\d+\s*(,\s*\d+\s*)+)';
        if ~isempty(regexpi(tmpstr,pat,'once'))
          error('Multidimensional waves not supported')
        end
      end

      % Check WAVES line for illegal characters (liberal wave name rules)
      if regexp(tmp,'[";:]+','once')
        error('Bad symbol on WAVES line')
      end
      numqts = numel(regexp(tmp,''''));

      % An even number of single quotes are expected on the WAVES line
      if ~isempty(numqts)
        if (numqts/2~=ceil(numqts/2))
          error('Bad symbol on WAVES line')
        end
      end

      % Get wave names
      count = 0;
      tmpstr = '';
      while (true)
        % Break from while loop if remaining WAVES line is empty
        if isempty(tmp)
          break
        end
        % Break from while loop if remaining WAVES line contains only delimiters
        if strcmp(regexprep(tmp,'['',\s]','','once'),'')
          break
        end
        % Find the position of the next standard and liberal wave names
        nxtstd = regexp(tmp,'[^'',\s]','once')-1;
        if isempty(nxtstd)
          nxtstd = inf;
        end
        nxtlib = regexp(tmp,'''','once')-1;
        if isempty(nxtlib)
          nxtlib = inf;
        end
        % Check for unrecognised characters trailing the WAVES line
        if isinf(nxtstd) && isinf(nxtlib)
          error('bad symbol on WAVES line')
        end
        % Check for illegal characters between wave names
        delim = tmp(1:min(nxtlib,nxtstd));
        delim = regexprep(delim,',',' ','once');
        if ~isempty(regexp(delim,'\S'))
          error('bad symbol on WAVES line')
        end
        % Remove wave name delimiter
        tmp(1:min(nxtlib,nxtstd)) = '';
        % Fetch next wavename
        if nxtstd < nxtlib
          % Standard wave names separated by a comma, tabs or spaces
          [tmpstr,tmp] = strtok(tmp,[',' char(9) char(32)]);
          % Check that the wave name contains no illegal characters
          if ~isempty(regexp(tmpstr,'\W'))
            error('bad symbol on WAVES line')
          end
          % Check that the wave name starts with a letter
          if isempty(regexpi(tmpstr(1),'[A-Z]'))
            error('bad symbol on WAVES line')
          end
        else
          % Liberal wave name rules implied by opening quotation mark
          [tmpstr,tmp] = strtok(tmp,'''');
          % Remove closing quotation mark left on WAVES line
          tmp(1) = '';
        end
        if isempty(names)
          names{1,1} = tmpstr;
        else
          names{end+1,1} = tmpstr;
        end
        % Step-up wave counter for current block of data
        count = 1+count;
      end

      % Check that the next line is the BEGIN keyword
      tmp = fgetl(fid);
      tmp = strjust(tmp,'left');
      if regexpi(tmp,'BEGIN','once') ~= 1
        error('expected BEGIN keyword after WAVES in Igor text file')
      end
      tmp(1:5) = [];
      if ~isempty(deblank(tmp))
        warning('Data detected on BEGIN line will be ignored')
      end

      % Get the data for the current block
      tmpdata = fscanf(fid,'%g',[count,inf])';

      % Add data block to the data array and get it's dimensions
      waves{block} = tmpdata;

      % Check that the next line is the END statement of a data block
      tmp=fgetl(fid);
      tmp = strjust(tmp,'left');
      if regexpi(tmp,'END','once') ~= 1
        error('expected END keyword after data block in Igor text file')
      end

    elseif ~isempty(regexpi(tmp,'X\s+SetScale','once'))
      tmp = regexprep(tmp,'X\s*','','once','ignorecase');
      n = numel(regexpi(tmp,'SetScale'));

      %% READ SETSCALE COMMANDS
      for i = 1:n
        [tmpstr,tmp] = strtok(tmp,';');
        tmparr = cell(0);
        for j = 1:5
          [tmparr{j},tmpstr] = strtok(tmpstr,[',' char(9) char(32)]);
        end
        [tmparr{j+1},tmpstr] = strtok(tmpstr,';');

        % SetScale for the X dimension
        if strcmpi(tmparr{2},'x')
          if isempty(delta) && ~isempty(startidx)
            error('Conflicting calls to SetScale for the x dimension')
          end
          % Assume per point scaling (/P) for the x dimension
          if isempty(regexpi(tmparr{1},'SetScale\s*/\s*P','once'))
            error('Per point scaling only is supported for the x dimension')
          end
          % num1 is starting index value
          num1 = tmparr{3};
          if isempty(startidx)
            startidx = eval(char(num1));
          else
            if startidx ~= eval(char(num1))
              error(['The starting index value for the X dimension is'...
                     ' not the same for all waves'])
            end
          end
          % num2 is the delta value
          num2 = tmparr{4};
          if isempty(delta)
            % Set the delta variable since it is undefined
            if strcmp(num2,'0')
              % Set delta to 1 if num2 is 0
              delta = 1;
            else
              delta = eval(char(num2));
            end
          else
            % Check that the value of num2 matches the delta variable
            if delta ~= eval(char(num2))
              error(['The delta value for the X dimension is not the'...
                     ' same for all waves'])
            end
          end
          % Read units of the X dimension
          if ~isempty(regexpi(tmparr{5},'".*"','once'))
            if isempty(xunit)
              % Set the X dimension unit since it is undefined
              xunit = tmparr{5}(2:end-1);
              if isempty(xunit)
                xunit = 's';
              end
            else
              if ~strcmp(tmparr{5}(2:end-1),xunit)
                error(['The unit for the X dimension is not the same'...
                       ' for all waves'])
              end
            end
          end

        % SetScale for the Y dimension
        elseif strcmpi(tmparr{2},'y')
          if ~isempty(regexpi(tmparr{1},'SetScale\s*/\s*P','once'))
            error('Inclusive or data full scaling only is supported for the Y dimension')
          end
          % Read units of the Y dimension
          if ~isempty(regexpi(tmparr{5},'".*"','once'))
            if isempty(xunit)
              % Set the X dimension unit since it is undefined
              xunit = tmparr{5}(2:end-1);
              startidx = waves{1}(1,1);
              delta = [];
              if isempty(regexp(tmparr{6},names{1},'once'))
                error('The first wave must define the X dimension units')
              end
              fprintf('ephysIO: The SetScale command did not define the X dimension.\n');
              fprintf('ephySIO: The first wave was assigned to the X dimension instead.\n');
            elseif isempty(yunit)
              % Set the Y dimension unit since it is undefined
              yunit = tmparr{5}(2:end-1);
              if isempty(yunit)
                yunit = 'au';
              end
            elseif ~strcmp(tmparr{5}(2:end-1),yunit)
              error(['The unit for the Y dimension is not the same for'...
                     ' all waves'])
            end
            if strcmp(xunit,'dat') || strcmp(yunit,'dat')
              error('Date/time waves are not supported')
            end
          end

        elseif strcmpi(tmparr{2},'t')
          error('SetScale for the chunks dimension is not supported')

        elseif strcmpi(tmparr{2},'z')
          error('SetScale for the layers dimension is not supported')
        end
      end

    else

      % Record Igor commands etc into notes cell array else skip the line
      if ~isempty(regexpi(tmp,'X\s+','once'))
        if isempty(notes)
          notes{1,1} = sprintf(tmp);
        else
          notes{end+1,1} = sprintf(tmp);
        end
      end

    end
  end
  fclose(fid);

  % Concatenate waves
  try
    waves = cell2mat(waves);
  catch
    error('The number of samples is not the same for all waves')
  end

  % Evaluate the X dimension
  if ~isempty(startidx) && ~isempty(delta)
    % Create x axis and add as first column of data array using startidx and delta
    waves = cat(2,[startidx:delta:startidx+(delta*(size(waves,1)-1))]',waves);
    xdiff = delta;
    if strcmp(xunit,'s')
      xdim{1} = 'Time';
    else
      xdim{1} = 'XWave';
    end
    names = cat(1,xdim,names);
  else
    dx = diff(waves(:,1));
    if any(diff(dx) > 1.192093e-07)
      xdiff = 0;
    else
      xdiff = dx(1);
    end
  end
  names = char(names);
  if isempty(xunit)
    xunit = '';
  end
  if isempty(yunit)
    yunit = '';
  end

  %--------------------------------------------------------------------------

  function ATFwrite (filename,array,xunit,yunit,names,notes)

  % Get number of columns from data array
  ncols = size(array,2);
  nopth = numel(notes);

  % Allocate unit and scale data appropriately
  if strcmp(xunit,'s')
    array(:,1) =  array(:,1)*1e+3;
    xunit = 'ms';
    names{1} = 'Time';
  elseif strcmp(xunit,'A')
    array(:,1) =  array(:,1)*1e+12;
    xunit = 'pA';
  elseif strcmp(xunit,'V')
    array(:,1) =  array(:,1)*1e+3;
    xunit = 'mV';
  else
    % Do nothing
  end
  if strcmp(yunit,'s')
    array(:,2:end) = array(:,2:end)*1e+3;
    yunit = 'ms';
  elseif strcmp(yunit,'A')
    array(:,2:end) = array(:,2:end)*1e+12;
    yunit = 'pA';
  elseif strcmp(yunit,'V')
    array(:,2:end) = array(:,2:end)*1e+3;
    yunit = 'mV';
  else
    % Do nothing
  end

  % Open a file for writing
  fid = fopen(filename,'w');

  % Print first header
  fprintf(fid,'ATF\t1.0\n');

  % Print second header
  fprintf(fid,'%d\t%d\n',nopth,ncols);

  % Print notes into optional header
  if nopth > 0
    for i = 1:nopth
      fprintf(fid,char(strcat(notes{i},'\n')));
    end
  end

  % Print column titles
  if ischar(names)
    names = cellstr(names);
  end
  fprintf(fid,'%s (%s)',deblank(names{1}),xunit);
  for i = 2:ncols
    fprintf(fid,'\t%s (%s)',deblank(names{i}),yunit);
  end
  fprintf(fid,'\n');

  % Print data matrix
  % X dimension: 64-bit double precision variables represent data to about
  % 15 significant figures (14 decimal places)
  dataformat = '%.14g';
  for i = 1:ncols-1
    % Y dimension: 32-bit single precision variables represent data to about
    % 7 significant figures (6 decimal places)
    dataformat = strcat(dataformat,'\t%.7g');
  end
  dataformat = strcat(dataformat,'\n');
  fprintf(fid,dataformat,array');

  % Close file identifier
  fclose(fid);

  %--------------------------------------------------------------------------

  function [array,xdiff,xunit,yunit,names,notes] = ATFread (filename)

  % Open a file for reading
  fid = fopen(filename,'r');

  % Read ATF keyword header
  tmp = fgetl(fid);
  tmpstr = strtok(tmp,char(9));
  if ~strcmp(tmpstr,'ATF')
    error('The file is not in tab-delimited Axon text format')
  end

  % Read the no. of optional headers and data columns from the 2nd header line
  tmp = fgetl(fid);
  tmpvec = sscanf(tmp,'%d');
  nopth = tmpvec(1);
  ncols = tmpvec(2);

  % Read all optional headers into notes array
  notes = cell(0);
  if nopth > 0
    for i = 1:nopth
      notes{i,1} = sprintf(fgetl(fid));
    end
  end
  notes = char(notes);

  % Get column titles and units
  tmp = fgets(fid);
  count = 0;
  tmpstr = '';
  names = cell(0);
  yunit = cell(0);
  while (true)
    % Break from while loop if remaining column titles line is empty
    if isempty(tmp)
      break
    end
    % Break from while loop if remaining column titles line contains
    % only trailing field separators, spaces and/or closing brackets
    if strcmp(regexprep(tmp,'[,\s)]',''),'')
      break
    end
    % Find the position of the next word character and next double quote
    nxtwrd = regexp(tmp,'[\w]','once')-1;
    if isempty(nxtwrd)
      nxtwrd = inf;
    end
    nxtqte = regexp(tmp,'["]','once');
    if isempty(nxtqte)
      nxtqte = inf;
    end
    % Remove field separator(s)
    tmp(1:min(nxtqte,nxtwrd)) = '';
    % Fetch next column title
    if nxtwrd < nxtqte
      % Find the next nearest field delimiter
      nxtdlm = regexp(tmp,'[),\t"\n]','once');
      % Unquoted column title (with units if applicable)
      tmpstr = deblank(tmp(1:nxtdlm-1));
      tmp(1:nxtdlm-1)='';
    else
      % Find the next nearest field delimiter
      nxtdlm = regexp(tmp,'[)"\n]','once');
      % Quoted column title (with units if applicable)
      tmpstr = deblank(tmp(1:nxtdlm-1));
      tmp(1:nxtdlm-1)='';
      % Remove closing quotation mark
      nxtqte = regexp(tmp,'["]','once');
      tmp(1:nxtqte)='';
    end
    if count == 0
      [tmpstr,rem] = strtok(tmpstr,'(');
      if isempty(rem)
        xunit = '';
      else
        xunit = regexprep(rem,'(\s*','');
      end
    else
      [tmpstr,rem] = strtok(tmpstr,'(');
      if isempty(rem)
        yunit{count,1} = '';
      else
        yunit{count,1} = regexprep(rem,'(\s*','');
      end
    end
    names{count+1,1} = deblank(strjust(tmpstr,'left'));
    % Step-up column counter
    count = 1+count;
  end
  if all(strcmp(yunit,yunit{1}))
    yunit = yunit{1};
  else
    error(['The unit for the Y dimension is not the same for'...
           ' all y columns'])
  end
  if numel(names) ~= ncols
    error('The number of column titles does not match the data array')
  end
  names = char(names);

  % Get the data array
  array = fscanf(fid,'%g',[ncols,inf])';

  % Evaluate X dimension
  dx = diff(array(:,1));
  if any(diff(dx) > 1.192093e-07)
    xdiff = 0;
  else
    xdiff = dx(1);
  end

%--------------------------------------------------------------------------

function [array,xdiff,xunit,yunit,names,notes] = MAload (filename,ch)

  % Read Metadata (requires readMeta.m file from 'acq4/pyqtgraph/metaarray.readMeta.m')
  metadata = readMeta(filename);
  if numel(metadata)<3
    error('The file is not a patch clamp recording trace')
  end

  % Read units for x and y axes
  xunit = metadata{2}.units(2);
  ch = ch+1;
  numChannels = numel(metadata{1}.cols);
  if ch > numChannels
    error('Channel number out of range');
  end
  yunit = metadata{1}.cols{ch}.units(2:end-1);

  % Read x-axis name into names array
  names = cell(0);
  names{1,1} = metadata{2}.name;

  % Read data
  % 0) Command channel
  % 1) Primary recording channel
  % 2) Secondary recording channel
  fprintf('Number of recording channels: %d\n',numChannels);
  for i = 1:numChannels
      fprintf('%d) %s\n',i-1,metadata{1}.cols{i}.name)
  end
  fprintf('loading channel %d...\n',ch-1)
  filepath = pwd;
  data = hdf5read(filename,'/data');
  l = size(data,1);
  array = zeros(l,2);
  array(:,1) = metadata{2}.values;
  array(:,2) = data(:,ch);
  names{2,1} = strcat('YWave',filepath(end-2:end));
  if filepath(end-3:end) == '/000'
    cd ..
    count = 1;
    exitflag = 0;
    while exitflag < 1
      data=[];
      dirname = strcat('00',num2str(count));
      dirname = dirname(end-2:end);
      if exist(dirname,'dir')
        cd(dirname);
        data = hdf5read(filename,'/data');
        array = cat(2,array,data(:,ch));
        names{count+2,1} = strcat('YWave',dirname);
        count = count + 1;
        cd ..
      else
        exitflag = 1;
      end
    end
  else
    % Do nothing
  end

  % Evaluate the X dimension
  dx = diff(array(:,1));
  if any(diff(dx) > 1.192093e-07)
    xdiff = 0;
  else
    xdiff = dx(1);
  end

  % Parse metadata into notes array
  notes = cell(0);
  obj = {'ClampState';'ClampParams';'DAQ';'Protocol'};
  for i=1:4
    if i == 2
      notes = cat(1,notes,[obj{i-1},'.',obj{i}]);
      key = fieldnames(metadata{3}.(obj{i-1}).(obj{i}));
      val = struct2cell(metadata{3}.(obj{i-1}).(obj{i}));
    elseif i == 3
      key = fieldnames(metadata{3}.(obj{i}));
      notes = cat(1,notes,[obj{3},'.',char(key(ch))]);
      tmp = struct2cell(metadata{3}.(obj{i}));
      key = fieldnames(tmp{ch});
      val = struct2cell(tmp{ch});
    else
      notes = cat(1,notes,[obj{i}]);
      key = fieldnames(metadata{3}.(obj{i}));
      val = struct2cell(metadata{3}.(obj{i}));
    end
    for j=1:numel(key)
      if ischar(val{j})
        notes = cat(1,notes,sprintf(['  ',char(key(j)),': ',val{j}]));
      elseif isnumeric(val{j})
        notes = cat(1,notes,sprintf(['  ',char(key(j)),': ',num2str(val{j})]));
      else
        % Do nothing
      end
    end
  end

  % Convert list of wave names to character array
  names = char(names);

%--------------------------------------------------------------------------

function [array,xdiff,xunit,yunit,names,notes] = H5load (filename,ch)

% Fetch recording attributes
rec_attrs = hdf5read(filename,'/description');
n_channels = rec_attrs.Data{1};
if ch > n_channels
  error('Channel number out of range');
end
info = hdf5info(filename);
groups = cellstr(char(info.GroupHierarchy.Groups.Name));
fprintf('Number of recording channels: %d\n',n_channels);
for i = 1:n_channels
    fprintf('%d) %s\n',i,groups{i}(2:end))
end
fprintf('loading channel %d...\n',ch)

% Fetch channel attributes
ch_attrs = hdf5read(filename,strcat(groups{ch},'/description'));
n_sections = ch_attrs.Data{1};

% Fetch data
names = cell(0);
charStr = num2str(zeros(numel(char(num2str(n_sections))),1))';
numChar = numel(charStr);
for i = 1:n_sections
   sectionNum = strcat(charStr,num2str(i-1));
   sectionNum = sectionNum(end-(numChar-1):end);
   path2sec = strcat(groups{ch},'/section_',sectionNum,'/data');
   data = hdf5read(filename,path2sec);
   if i==1
     nrows = numel(data);
     array = zeros(nrows,n_sections);
     array(:,1) = data;
     path2attrs = strcat(groups{ch},'/section_',sectionNum,'/description');
     sec_attrs = hdf5read(filename,path2attrs);
     % Parse section attributes
     dt = sec_attrs.Data{1};
     xdiff = dt;
     xunit = sec_attrs.Data{2}.Data;
     yunit = sec_attrs.Data{3}.Data;
   else
     % Check that the section attributes conform to the first section
     if numel(data)~=nrows
       error('The number of samples is not the same for all sections')
     end
     path2attrs = strcat(groups{ch},'/section_',sectionNum,'/description');
     sec_attrs = hdf5read(filename,path2attrs);
     if sec_attrs.Data{1}~=dt
       error(['The delta value for the X dimension is not the'...
              ' same for all sections'])
     end
     if sec_attrs.Data{2}.Data~=xunit
       error(['The unit for the X dimension is not the same for'...
              ' all sections'])
     end
     if sec_attrs.Data{3}.Data~=xunit
       error(['The unit for the Y dimension is not the same for'...
              ' all sections'])
     end
     array(:,i) = data;
   end
   names{i+1,1} = strcat('Section_',sectionNum);
   data = [];
end

% Construct X dimension and add it to the data array
x = dt*[0:nrows-1]';
array = cat(2,x,array);
names{1} = 'Time';

% Scale dimenions
if xunit=='ms'
  array(:,1) = 1e-3*array(:,1);
  xunit = 's';
end
if yunit=='pA'
  array(:,2:end) = 1e-12*array(:,2:end);
  yunit = 'A';
elseif yunit=='mV'
  array(:,2:end) = 1e-3*array(:,2:end);
  yunit = 'V';
end

% Create empty notes cell array
notes = cell(0);

%--------------------------------------------------------------------------

function [array,xdiff,xunit,yunit,names,notes] = IBWload (filename)

  % IBWread helper function to get data from Igor binary wave
  S = IBWread (filename);

  % Create x axis and add as first column of data array
  array = cat(2,[S.x0:S.dx:S.x0+(S.dx*(S.Nsam-1))]',S.y);
  xdiff = S.dx;

  % Create blank notes variable
  try
    notes = S.WaveNotes;
    tmparr = strread(S.WaveNotes,'%s','delimiter','\n:');
  catch
    notes = {};
    tmparr = {};
  end

  % Get X dimension units
  try
    xunit = S.waveHeader.xUnits;
    if isempty(xunit)
      error('invoke catch statement')
    end
  catch
    try
      xunit = tmparr(find(strcmp(tmparr,'xLabel'))+1);
    catch
      xunit = 's';
    end
  end

  % Get Y dimension units
  try
    yunit = S.waveHeader.dataUnits;
    if isempty(yunit)
      error('invoke catch statement')
    end
  catch
   try
     yunit = tmparr(find(strcmp(tmparr,'yLabel'))+1);
   catch
     yunit = 'au';
   end
  end

  % Get X and Y dimension names
  names = cell(2,1);
  if strcmp(xunit,'s')
    names{1} = 'Time';
  else
    names{1} = 'XWave';
  end
  names{2} = S.bname;

%--------------------------------------------------------------------------

function [array,xdiff,xunit,yunit,names,notes] = ABF2load (filename,ch)

  [d,si,h] = abfload(filename);

  % Evaluate the recording channels
  fprintf('Number of recording channels: %d\n',h.nADCNumChannels);
  for i = 1:h.nADCNumChannels
    fprintf('%d) %s\n',i,h.recChNames{i});
  end
  fprintf('loading channel %d...\n',ch)
  if ch > h.nADCNumChannels
    error('Channel number out of range')
  end

  % Get data for the selected recording channel
  if h.nOperationMode == 3
    array = d(:,ch);
  else
    array = reshape(d(:,ch,:),h.dataPtsPerChan,h.lActualEpisodes);
  end
  si = si*1e-6;
  array = cat(2,[0:si:si*(h.dataPtsPerChan-1)]',array);
  xdiff = si;
  xunit = 's';
  yunit = strjust(h.recChUnits{ch},'left');
  notes = '';

  % Assign question marks to character array of column names
  ncols = size(array,2);
  names = char(zeros(ncols,1)+63);
