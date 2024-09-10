from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain.agents.format_scratchpad.openai_tools import format_to_openai_tool_messages
from langchain.agents.output_parsers.openai_tools import OpenAIToolsAgentOutputParser
from langchain_core.messages import AIMessage, HumanMessage
from langchain.agents import AgentExecutor

class Agent(object):

    def __init__(self, llm_provider, llm_name, api_key, temperature, tools, sys_prompt, initial_message) -> None:
        """
        Initializes the Agent object.

        Parameters:
        - llm_provider: The language model provider (e.g., 'openai').
        - llm_name: The name of the language model.
        - api_key: The API key for the language model provider.
        - temperature: The temperature parameter for generating responses.
        - tools: The tools to bind to the language model.
        - sys_prompt: The system prompt for the conversation.
        - initial_message: The initial message in the conversation.

        Raises:
        - ValueError: If the LLM provider is not supported.
        """

        # Set the instance variables
        self.provider = llm_provider
        self.llm_name = llm_name.lower()
        self.temperature = temperature
        self.tools = tools
        self.sys_prompt = sys_prompt
        self.chat_history = []
        MEMORY_KEY = "chat_history"

        # Add the initial message to the chat history
        self.chat_history.extend([AIMessage(content=initial_message)])

        # Check if the language model provider is OpenAI
        if self.provider == 'openai':   
            # Initialize the ChatOpenAI object with the specified parameters
            self.llm = ChatOpenAI(model=self.llm_name, openai_api_key=api_key, temperature=self.temperature)

        else:
            # Raise an error if the language model provider is not supported
            raise ValueError('LLM provider not supported')
        
        # Bind the tools to the language model
        self.llm_with_tools = self.llm.bind_tools(self.tools)

        # Create the prompt template for the conversation
        self.prompt = ChatPromptTemplate.from_messages(
            [
            ('system', sys_prompt),
            MessagesPlaceholder(variable_name=MEMORY_KEY),
            ("user", "{input}"),
            MessagesPlaceholder(variable_name="agent_scratchpad"),
            ]
        )

        # Check if the language model provider is OpenAI
        if self.provider == 'openai':
            # Define the agent pipeline using the prompt template, language model with tools, and output parser
            self.agent = (
            {
                "input": lambda x: x["input"],
                "agent_scratchpad": lambda x: format_to_openai_tool_messages(x["intermediate_steps"]),
                "chat_history": lambda x: x["chat_history"],
            }
            | self.prompt
            | self.llm_with_tools
            | OpenAIToolsAgentOutputParser()
            )

        # Create the agent executor with the agent pipeline and tools
        self.agent_executor = AgentExecutor(agent=self.agent, tools=tools, verbose=True)


    def chat(self, message: str) -> str:
        """
        Sends a message to the agent and returns the agent's response.

        Parameters:
        - message: The message to send to the agent.

        Returns:
        - The agent's response as a string.
        """

        result = self.agent_executor.invoke({"input": message, "chat_history": self.chat_history})

        self.chat_history.extend([HumanMessage(content=message), AIMessage(content=result["output"])])

        return result["output"]
