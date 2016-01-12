%% @author Jörg Brünecke <devel@bloerg.de>
%% @doc FITS format input/output routines
%% The module provides functions for reading, writing and parsing 
%% FITS files
%% @reference <a href="http://fits.gsfc.nasa.gov/fits_documentation.html">NASA FITS Documentation</a>

-module(fitserl).
-export([load_fits_file/1]).
-export([get_hdus/1, get_hdu_header/1, get_hdu_header/2]).
-export([get_hdu_data/1, get_hdu_data/2, get_hdu_data/3, find_hdu_data/1]).
-export([parse_binary_table/2, binary_table_extract_rows/2, binary_table_extract_fields/2]).


%% @doc Load a FITS file and return the binary content
load_fits_file(File_path) 
when 
    is_list(File_path) ->

    case file:read_file(File_path) of 
        {error, enoent} -> {error, file_does_not_exist};
        {error, eaccess} -> {error, access_to_file_permitted};
        {error, eisdir} -> {error, is_not_a_file};
        {error, enotdir} -> {error, path_error};
        {error, enomem} -> {error, not_enough_memory};
        {ok, Binary_data} ->
            {ok, Binary_data}
    end.

%% @doc Tries to find out the type of a value encoded as
%% a string. If no type can be determined the input string is returned
try_type_cast(Input_string) 
when is_list(Input_string) ->
    Value_string = string:strip(Input_string),
    case hd(Value_string) of
        %% a string starting with single quote (ASCII 39)
        39 -> 
            % try to get the string within the single quotes without blanks
            Temp_string = string:substr(Value_string, string:chr(Value_string, 39)+1),
            string:strip(
                string:substr(
                    Temp_string,
                    1,
                    case string:chr(Temp_string, 39) -1 of 
                        -1 -> length(Temp_string);
                        Else -> Else
                    end
                )
            );
        %% possibly a number
        _ ->
            case string:str(Value_string, ".") of
                % no "." can be found -> possibly integer
                0 -> 
                    try string:to_integer(Value_string) of
                        {Value,[]} -> Value;
                        %String contains characters after the number
                        {_Value, _Rest} -> Value_string
                    catch {error, no_integer} ->
                        Value_string
                    end;
                % a "." can be found -> possibly float
                _Else -> 
                    try string:to_float(Value_string) of
                        {Value, []} -> Value;
                        %String contains characters after the number
                        {_Value, _Rest} -> Value_string
                    catch {error, no_float} ->
                        Value_string
                    end
            end
    end.

%% @doc helper function for starting parse_header/2
%% Input is Binary_data containing an hdu starting at offset 0
get_hdu_header(Binary) when is_binary(Binary) ->
    get_hdu_header(Binary, []).

%% @doc helper function for starting parse_header/2
%% with the number of the hdu as input (starting at 1)
get_hdu_header(Hdu_number, Fits) 
when is_integer(Hdu_number),
     is_binary(Fits) ->
        % get the byte offset in Fits of the hdu with number Hdu_number
        case lists:keyfind(Hdu_number, 1, get_hdus(Fits)) of
            {Hdu_number, Pos} -> 
                Next_hdu_number = Hdu_number + 1,
                % get the length of the hdu with Hdu_number from the
                % difference of (the staring position of the next hdu or
                % the end position of Fits) and the start position of
                % the hdu
                Len =
                    case lists:keyfind(Next_hdu_number, 1, get_hdus(Fits)) of
                        {Next_hdu_number, Byte_end} -> Byte_end - Pos;
                        false -> byte_size(Fits) - Pos
                    end,
                % start parsing
                get_hdu_header(binary_part(Fits, Pos, Len));
            % there is no hdu with Hdu_number
            false -> {error, out_of_hdu_range}
        end;

%% @doc recursively parses a binary for Key, Value, Comment tuples
%% stops parsing on the occurence of a line starting with "END"
%% Returns a list of {Keyword, Value, Comment} tuples for the header
get_hdu_header(Binary, Plain_text_header) ->
    case Binary of 
        % split off the next/last line of 80 characters
        <<Line:80/binary, Rest/binary>> -> ok;
        <<Line:80/binary>> -> Rest = <<>>
    end,
    case binary_part(Line, 0, 3) of
        % end of hdu header, return result
        <<"END">> -> Plain_text_header;
        % go on processing keyword record
        _ -> 
            % comments are separated by "/"
            case binary:split(Line, <<"/">>) of
                % FIXME: What if the Value or comment contain one or
                % several "/"?
                [Key_and_value, Comment] ->
                    Comment_string = string:strip(binary_to_list(Comment));
                [_Else] -> 
                    Comment_string = [],
                    Key_and_value = Line
            end,
            % keywords and values are separated by "="
            case binary:split(Key_and_value, <<"=">>) of
                [Key_word, Value] ->
                    Key_word_string = string:strip(binary_to_list(Key_word)),
                    Value_string = string:strip(binary_to_list(Value));
                [String_without_equal_sign] ->
                    %% see http://fits.gsfc.nasa.gov/fits_primer.html:
                    %% Some keywords, (e.g., COMMENT and HISTORY) are 
                    %% not followed by an equals sign and in that case 
                    %% columns 9 - 80 of the record may contain any 
                    %% string of ASCII text.
                    <<Key_word:8/binary, Value/binary>> = String_without_equal_sign, 
                    Key_word_string = string:strip(binary_to_list(Key_word)),
                    Value_string = string:strip(binary_to_list(Value))
            end,
            get_hdu_header(Rest, [{Key_word_string, try_type_cast(Value_string), Comment_string}|Plain_text_header])
    end.


%% @doc helper function fo rstarting find_hdu_data/2
find_hdu_data(Binary) when is_binary(Binary)-> 
    find_hdu_data(Binary,0).

%% @doc recursively walks through binary trying to find a line beginning
%% with "END". Then returns the byte offset of 2880 bytes * N + 1 after
%% the "END". This is the next possible position for a data section in 
%% an hdu.
find_hdu_data(Binary, Blocks_of_80_bytes) ->
    case Binary of
        % split off the next/last line of 80 characters
        <<Line:80/binary, Rest/binary>> -> ok;
        <<Line:80/binary>> -> Rest = <<>>;
        _ -> 
            Line = <<>>,
            Rest = <<>>
    end,
    case binary_part(Line, 0, min(3, byte_size(Line))) of
        % end of hdu header, return result
        <<"END">> -> 
            Possible_data_start_pos = (Blocks_of_80_bytes * 80) 
             + (2880 - Blocks_of_80_bytes * 80 rem 2880),
             Possible_data_start_pos;
        <<>> -> {error, malformed_hdu};
        _Else -> 
            find_hdu_data(Rest, Blocks_of_80_bytes + 1)
    end.
    

%% @doc helper function for starting get_hdu_data/1
%% with the number of the hdu as input (starting at 1)
get_hdu_data(Hdu_number, Fits)
when is_integer(Hdu_number),
     is_binary(Fits) ->
        % get the byte offset in Fits of the hdu with number Hdu_number
        case lists:keyfind(Hdu_number, 1, get_hdus(Fits)) of
            {Hdu_number, Pos} -> 
                Next_hdu_number = Hdu_number + 1,
                % get the length of the hdu with Hdu_number from the
                % difference of (the staring position of the next hdu or
                % the end position of Fits) and the start position of
                % the hdu
                Len =
                    case lists:keyfind(Next_hdu_number, 1, get_hdus(Fits)) of
                        {Next_hdu_number, Byte_end} -> Byte_end - Pos;
                        false -> byte_size(Fits) - Pos
                    end,
                % start parsing
                get_hdu_data(binary_part(Fits, Pos, Len));
            % there is no hdu with Hdu_number
            false -> {error, out_of_hdu_range}
        end.
                
%% @doc helper function for starting get_hdu_data/3
%% Input is a Binary containing an hdu
%% calls find_hdu(Binary) to get the possible offset of the
%% first data section within the Binary
get_hdu_data(Binary) when is_binary(Binary)->
    get_hdu_data(Binary, find_hdu_data(Binary), <<>>).

%% @doc get the binary data out of an hdu
%% returns <<>> if there is no data and the binary data representation
%% if data exists
get_hdu_data(Binary, Start_offset, Data) ->
    case Binary of 
        % end of binary
        <<>> -> Data;
        % next hdu starts
        <<_:(Start_offset)/binary,"XTENSION",_/binary>> -> Data;
        % next one should not happen because SIMPLE /can/ only occur at
        % the beginning of a fits file
        <<_:(Start_offset)/binary,"SIMPLE",_/binary>> -> Data; 
        %within the data section
        <<_:(Start_offset)/binary,Block_of_2880_bytes:2880/binary, Rest/binary>> ->
            get_hdu_data(Rest, 0, <<Data/binary, Block_of_2880_bytes/binary>>);
        %hdu has no data section
        _Else -> <<>>
    end.


%% @doc recursively find all hdus and return a list of 
% {Hdu_number, Header_start_byte_offset} tuples
% returns [] if no hdu is found
get_hdus(Fits) ->
    get_hdus(Fits, 0, []).
get_hdus(Fits, Block_count, Result) ->
    case Fits of 
        <<Key_word:8/binary, _:2872/binary, Rest/binary>> -> ok;
        <<>> -> 
            Key_word = <<>>,
            Rest = <<>>;
        <<Rest/binary>> -> Key_word = <<>>  % this acually means a truncated hdu
    end,
    case Key_word of
        % end of Fits
        <<>> -> 
            lists:reverse(Result);
        % primary hdu
        <<"SIMPLE",_/binary>> ->
            New_result = [{length(Result) + 1, Block_count*2880}|Result],
            get_hdus(Rest, Block_count + 1, New_result);
        % additional hdu
        <<"XTENSION",_/binary>> ->
            New_result = [{length(Result) + 1, Block_count*2880}|Result],
            get_hdus(Rest, Block_count + 1, New_result);
        % somewhere in between
        _Else -> 
            get_hdus(Rest, Block_count + 1, Result)
    end.

%% @doc extract the Rows from a hdu data binary
%% The Row_length parameter should be taken from the NAXIS1 value of
%% the corresponding header
%% Returns [] or a List of binaries with each binary being a row
binary_table_extract_rows(Data, Row_length) ->
    binary_table_extract_rows(Data, Row_length, []).
binary_table_extract_rows(Data, Row_length, Result) ->
    case Data of
        <<>> -> lists:reverse(Result);
        <<Row:Row_length/binary, Rest/binary>> ->
            binary_table_extract_rows(Rest, Row_length, [Row|Result]);
        _Else -> {error, malformed_hdu_data}
    end.



binary_table_extract_fields(Row_data, TFORM)
when is_binary(Row_data), is_list(TFORM) ->
    binary_table_extract_fields(Row_data, TFORM, []).
binary_table_extract_fields(Row_data, TFORM, Result) ->
    Field_types =
        [
        {"L", {logical, 1,
             fun(Binary) ->
                case Binary of
                    <<"T">> -> true;
                    <<"F">> -> false;
                    _Else -> {error, value_not_allowed}
                end
            end}
        },
        {"X", {bit_array, 1,
             fun(Binary) -> Binary %FIXME: I don't understand the definition
             end
            }
        },
        {"B", {unsigned_8_bit_integer, 1,
             fun(Binary) -> 
                case Binary of 
                    <<Integer8bu:8/big-unsigned-integer-unit:1>> -> Integer8bu;
                    _Else -> {error, value_not_allowed}
                end
            end}
        },
        {"I", {signed_16_bit_integer, 2,
             fun(Binary) ->
                case Binary of 
                    <<Integer16bs:16/big-signed-integer>> -> Integer16bs;
                    _Else -> {error, value_not_allowed}
                end
            end}
        },
        {"J", {signed_32_bit_integer, 4,
             fun(Binary) ->
                case Binary of
                    <<Integer32bs:32/big-signed-integer>> -> Integer32bs;
                    _Else -> {error, value_not_allowed}
                end
            end}
        },
        {"A", {character_string, 1,
             fun(Binary) ->
                case Binary of 
                    <<0>> -> [];
                    <<_Character>> -> binary_to_list(Binary);
                    _Else -> {error, value_not_allowed}
                end
            end}
        },
        {"E", {single_precision_float, 4,
             fun(Binary) ->
                case Binary of
                    <<Float32bs:32/big-float>> -> Float32bs;
                    _Else -> {error, value_not_allowed}
                end
            end}
        },
        {"D", {double_precision_float, 8,
             fun(Binary) ->
                case Binary of
                    <<Float64bs:64/big-float>> -> Float64bs;
                    _Else -> {error, value_not_allowed}
                end
            end}
        },
        {"C", {single_precision_complex, 8,
             fun(Binary) ->
                case Binary of
                    <<Real32bs:32/big-float,Imaginary32bs:32/big-float>> -> {Real32bs, Imaginary32bs};
                    _Else -> {error, value_not_allowed}
                end
            end}
        },
        {"M", {double_precision_complex, 16,
             fun(Binary) ->
                case Binary of
                    <<Real64bs:64/big-float,Imaginary64bs:64/big-float>> -> {Real64bs, Imaginary64bs};
                    _Else -> {error, value_not_allowed}
                end
            end}
        },
        {"P", {array_descriptor, 8,
             fun(Binary) ->
                case Binary of
                    <<Integer32b_1:32/big-integer,Integer32b_2:32/big-integer>> -> {Integer32b_1, Integer32b_2};
                    _Else -> {error, value_not_allowed}
                end
            end}
        }
        ],
    case TFORM of 
        [] -> 
            Field_type = [],
            Byte_length = 0,
            Handle_field_fun = fun() -> [] end;
        _ ->
            {_,Field_type,_} = hd(TFORM),
            {_Key, {_Type, Byte_length, Handle_field_fun}} = lists:keyfind(Field_type, 1, Field_types)
    end,
    case Row_data of
        <<>> -> lists:reverse(Result);
        <<Field:Byte_length/binary, Rest/binary>> -> 
            binary_table_extract_fields(Rest,tl(TFORM), [Handle_field_fun(Field)|Result]);
        _Else -> {error, malformed_row_data}
    end.
    

parse_binary_table(Data, Header) ->

    {_, BITPIX, _}  = lists:keyfind("BITPIX", 1, Header),
    {_, NAXIS, _}   = lists:keyfind("NAXIS", 1, Header),
    {_, NAXIS1, _}  = lists:keyfind("NAXIS1", 1, Header),
    {_, NAXIS2, _}  = lists:keyfind("NAXIS2", 1, Header),
    {_, PCOUNT, _}  = lists:keyfind("PCOUNT", 1, Header),
    {_, GCOUNT, _}  = lists:keyfind("GCOUNT", 1, Header),
    {_, TFIELDS, _} = lists:keyfind("TFIELDS", 1, Header),
    TFORM = [ lists:keyfind("TFORM" ++ integer_to_list(N), 1, Header) 
              || N <- lists:seq(1, TFIELDS)
            ],
    TTYPE = [ lists:keyfind("TTYPE" ++ integer_to_list(N), 1, Header) 
              || N <- lists:seq(1, TFIELDS)
            ],
    TUNIT = [ lists:keyfind("TUNIT" ++ integer_to_list(N), 1, Header) 
              || N <- lists:seq(1, TFIELDS)
            ],
    TNULL = [ lists:keyfind("TNULL" ++ integer_to_list(N), 1, Header) 
              || N <- lists:seq(1, TFIELDS)
            ],
    TDISP = [ lists:keyfind("TDISP" ++ integer_to_list(N), 1, Header) 
              || N <- lists:seq(1, TFIELDS)
            ],
    TSCAL = [ lists:keyfind("TSCAL" ++ integer_to_list(N), 1, Header) 
              || N <- lists:seq(1, TFIELDS)
            ],
    TZERO = [ lists:keyfind("TSCAL" ++ integer_to_list(N), 1, Header) 
              || N <- lists:seq(1, TFIELDS)
            ]
    

.
    
    
