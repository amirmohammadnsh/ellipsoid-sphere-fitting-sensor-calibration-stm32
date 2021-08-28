/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * <h2><center>&copy; Copyright (c) 2021 STMicroelectronics.
  * All rights reserved.</center></h2>
  *
  * This software component is licensed by ST under BSD 3-Clause license,
  * the "License"; You may not use this file except in compliance with the
  * License. You may obtain a copy of the License at:
  *                        opensource.org/licenses/BSD-3-Clause
  *
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"
#include "gpio.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include "ellipsoidfit.h"
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

/* USER CODE BEGIN PV */
#define numberOfSamples 100
#define numberOfAxis 3



const float32_t samples[numberOfSamples*numberOfAxis] ={
		294,-1508,-2,
		299,-1518,1,
		299,-1518,1,
		300,-1523,-11,
		297,-1538,8,
		310,-1543,7,
		310,-1543,28,
		320,-1544,28,
		329,-1551,36,
		368,-1562,62,
		407,-1557,82,
		407,-1557,82,
		490,-1545,149,
		549,-1553,194,
		594,-1538,246,
		652,-1503,310,
		652,-1503,310,
		680,-1457,356,
		707,-1401,381,
		730,-1321,380,
		727,-1247,357,
		727,-1247,357,
		695,-1155,313,
		643,-1118,250,
		592,-1085,134,
		592,-1078,44,
		513,-1078,44,
		405,-1118,-46,
		283,-1201,-113,
		147,-1277,-142,
		147,-1277,-142,
		-28,-1357,-119,
		-136,-1444,-62,
		-249,-1504,33,
		-334,-1554,117,
		-334,-1554,117,
		-427,-1567,197,
		-462,-1577,284,
		-513,-1547,338,
		-538,-1516,339,
		-538,-1516,339,
		-549,-1503,345,
		-540,-1499,324,
		-511,-1516,280,
		-471,-1519,230,
		-471,-1519,230,
		-409,-1574,204,
		-316,-1625,164,
		-230,-1664,136,
		-230,-1664,136,
		-114,-1687,124,
		30,-1717,136,
		139,-1727,142,
		248,-1728,184,
		248,-1728,184,
		362,-1701,210,
		451,-1639,216,
		526,-1580,227,
		590,-1474,219,
		590,-1474,219,
		609,-1372,196,
		599,-1289,149,
		566,-1204,79,
		482,-1142,11,
		482,-1142,11,
		378,-1117,-59,
		265,-1122,-102,
		144,-1172,-143,
		144,-1172,-143,
		30,-1238,-152,
		-100,-1303,-125,
		-217,-1372,-82,
		-301,-1451,28,
		-301,-1451,28,
		-386,-1490,109,
		-456,-1501,180,
		-507,-1469,229,
		-523,-1454,245,
		-523,-1454,245,
		-555,-1423,256,
		-550,-1421,225,
		-508,-1432,184,
		-445,-1478,141,
		-445,-1478,141,
		-384,-1512,121,
		-304,-1569,74,
		-204,-1604,52,
		-204,-1604,52,
		-87,-1642,44,
		35,-1683,63,
		141,-1682,73,
		241,-1664,96,
		241,-1664,96,
		349,-1638,111,
		424,-1592,127,
		490,-1552,146,
		552,-1447,141,
		552,-1447,141,
		581,-1370,136,
		567,-1266,99};
double offset[3]={0};
double gain[3] = {0};
double rotM[9] = {0};
double *results;
unsigned long start =0 ;
unsigned long end = 0;
unsigned long elapsed = 0;

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
/* USER CODE BEGIN PFP */

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{
  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  /* USER CODE BEGIN 2 */
  CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
  DWT->CYCCNT = 0;
  DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;
  start = DWT -> CYCCNT;

  results = ellipsoid_fit(&samples,numberOfSamples*numberOfAxis);
  for(int i = 0 ; i < 3; i++)
	  offset[i] = *(results + i);

  for(int i = 3; i <6; i++)
	  gain[i-3] =*(results + i);

  for(int i = 6; i <15; i++)
	  rotM[i-6] =*(results + i);

  refine_3D_fit(&gain,&rotM);
  end = DWT -> CYCCNT;
  elapsed = end -start;
  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  while (1)
  {
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSE;
  RCC_OscInitStruct.HSEState = RCC_HSE_ON;
  RCC_OscInitStruct.HSEPredivValue = RCC_HSE_PREDIV_DIV1;
  RCC_OscInitStruct.HSIState = RCC_HSI_ON;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSE;
  RCC_OscInitStruct.PLL.PLLMUL = RCC_PLL_MUL9;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }
  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV2;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_2) != HAL_OK)
  {
    Error_Handler();
  }
}

/* USER CODE BEGIN 4 */

/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */

/************************ (C) COPYRIGHT STMicroelectronics *****END OF FILE****/
