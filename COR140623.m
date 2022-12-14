strD= 'd:\Work\MATLAB\Electronic\SPLIT002.txt';
strC= '\\Tamara\D\Work\ATS\temba2013\1406TEST\Termochain\N3\Test\coef#3_ 2cor_140636.txt';

SLoad= struct('strDir',strD, 'col',struct( 'FirstRow',2, ...
  'inDate','', 'DateForm',{'yyyy m d H M S'}, ...
  'header'       , '*yyyy mm dd HH MM SS	P0_dBar T00 0 Power_V P1_dBar P2_dBar T%02d', ...
  'headerReplace', ['yyyy mm dd HH MM SS  Time P0_dBar P1_dBar P2_dBar Power_V T00 T%02d', ...
  char(13), 'Date1 Date2 Date3 Date4 Date5 Date6 Time P0_dBar P1_dBar P2_dBar Power_V T00 T%02d'], ...
  'ColSeparator',char(9), 'notWarningLostFields',{'Time'}));

DATA= termochain(SLoad,strC);