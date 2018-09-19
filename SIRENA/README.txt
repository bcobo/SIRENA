1. Simular con xifusim
2. Reconstruir con el SIXTE de 9May o con el nuevo pero cambiando el fichero 'testriggerfile.c'(*)

(*) SIXTE-9May: Funciona con ADC variable                        |
    SIXTE-nuevo(18/09/2018): Funcionan con la columna ADC fija   | =>  Por eso para reconstruir hay que usar el SIXTE de 9May o el nuevo cambiando 'testriggerfile.c'
    xifusim: Crea la columna ADC de longitud variable            |



'testriggerfile.c' (9May:OK, Nuevo(18/09/2018):Mal):
   ADC de longitud variable: 
	LONGLONG rec_trigsize;
    	LONGLONG offset;
	fits_read_descriptll(file->fptr,file->trigCol,file->row,&rec_trigsize,&offset,status);
   ADC de longitud fija:
	LOLONGLONG rec_trigsize;
    	LONGLONG col_width;
	int adc_col_typecode;
    	fits_get_coltypell(file->fptr,file->trigCol,&adc_col_typecode,
		       &rec_trigsize,&col_width,status);

   Si usa 'fits_get_coltypell' para leer una columna de longitud variable, rec_trigsize=1 => record->trigger_size=1 en integraSIRENA.cpp => reconstruct_init->pulse_length > record->trigger_size => ERROR="Pulse length is larger than record size"



'tes_simulation.c' (9May:OK):
   Había que comentar o no la línea donde se define el OVERSAMPLING para poder hacer jitter o no
   No lo vamos a usar ahora porque usaremos XIFUSIM en lugar de TESSIM para simular



'testconstpileup.c' (Nuevo(18/09/2018):OK):
   Se puede seguir usando la opción offset=-1



'tessim_trigger.c' (Nuevo(18/09/2018):OK):
   Hubo que cambiar 'data->stream->time=time - data->preBufferSize*tes->delta_t;' por 'data->stream->time=time-data->preBufferSize*tes->delta_t*tes->decimate_factor;' (Ver correos con PP y Christian)

