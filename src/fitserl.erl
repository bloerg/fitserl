-module(fitserl).
%% @author Jörg Brünecke <devel@bloerg.de>
%% @doc FITS format input/output routines
%% The module provides functions for reading, writing and parsing 
%% FITS files
%% @reference <a href="http://fits.gsfc.nasa.gov/fits_documentation.html">NASA FITS Documentation</a>

-export([load_fits_file/1]).
-export([get_hdus/1, parse_hdu/1]).



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
    Value_binary = list_to_binary(Value_string),
    case [binary:first(Value_binary), binary:last(Value_binary)] of
        %% a string in single quotes
        [<<"'">>,<<"'">>] -> Value_string;
        %% possibly a number
        [_,_] ->
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
parse_hdu(Binary) ->
    parse_hdu(Binary, []).
%% @doc Parses a binary for Key, Value, Comment tuples
%% stops parsing on the occurence of a line starting with "END"
%% Returns a List of {Keyword, Value, Comment} tuples for the header
parse_hdu(Binary, Plain_text_header) ->
    <<Line:80/binary, Rest/binary>> = Binary,
    case binary_part(Line, 0, 3) of
        <<"END">> -> Plain_text_header;
        _ -> 
            case binary:split(Line, <<"/">>) of
                % FIXME: What if the Value or comment contain one or
                % several "/"?
                [Key_and_value, Comment] ->
                    Comment_string = string:strip(binary_to_list(Comment));
                [_Else] -> 
                    Comment_string = [],
                    Key_and_value = Line
            end,
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
            parse_hdu(Rest, [{Key_word_string, try_type_cast(Value_string), Comment_string}|Plain_text_header])
    end.


%% @doc find all hdus and return a list of 
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
            New_result = [{Block_count+1, Block_count*2880}|Result],
            get_hdus(Rest, Block_count + 1, New_result);
        % additional hdu
        <<"XTENSION",_/binary>> ->
            New_result = [{Block_count+1, Block_count*2880}|Result],
            get_hdus(Rest, Block_count + 1, New_result);
        % somewhere in between
        _Else -> 
            get_hdus(Rest, Block_count + 1, Result)
    end.

